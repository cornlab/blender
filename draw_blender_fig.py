import svgwrite
import argparse
import csv

# Stacia Wyman 21 July 2019
# Python script to draw alignments of BLENDER hit to the target sequence 
# Adapted from visualize.py from https://github.com/aryeelab/guideseq


boxWidth = 10
box_size = 16

colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', 'R': '#FFA500'}

REQUIRED_COLUMNS = (
    'Chr:Start-End', 'Discoscore', 'PAM', 'Guide sequence', 'Mismatches'
)
OPTIONAL_METRICS = ('DPM', 'normDPM')


class SiteList(list):
    """Parsed sites plus the optional metric columns declared by the input."""

    def __init__(self, available_metrics=()):
        super().__init__()
        self.available_metrics = tuple(available_metrics)

def parse_arguments():
    parser = argparse.ArgumentParser(prog="blender figure", description = "Make pretty figures from blender output")
    parser.add_argument('-i', '--infile', type=str, default=None, required=True, help='blender output to draw')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output filename')
    parser.add_argument('-g', '--guide', type=str, required=True, help='Guide RNA sequence (5\' to 3\')')
    parser.add_argument('-p', '--pam', metavar='PAM', required=True, help='PAM to output (5\' to 3\')')
    parser.add_argument('-l', '--location', type=int, default=3, help="PAM location (5\' or 3\' of the guide RNA). Default 3\'")
    parser.add_argument('-t', '--title', type=str, help="Optional title to add to the plot")
    parser.add_argument('--notext', action='store_true')
    args = parser.parse_args()
    return args

# Read in filtered off-target list, grab relevant columns, sort by disco score, then mismatches, 
# then genomic coordinate
def parseSitesFile(infile):
    with open(infile, 'r', encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames is None:
            raise ValueError('input is empty or has no header')
        missing = [column for column in REQUIRED_COLUMNS
                   if column not in reader.fieldnames]
        if missing:
            raise ValueError(
                'input is missing required columns: ' + ', '.join(missing)
            )
        available_metrics = [metric for metric in OPTIONAL_METRICS
                             if metric in reader.fieldnames]
        offtargets = SiteList(available_metrics)

        for line_number, row in enumerate(reader, start=2):
            if not any(isinstance(value, str) and value.strip()
                       for value in row.values()):
                continue
            try:
                disco_score = int(row['Discoscore'].strip())
                mismatches = int(row['Mismatches'].strip())
                offtarget_sequence = row['Guide sequence'].strip()
                pam = row['PAM'].strip()
                coords = row['Chr:Start-End'].strip()
                metric_values = {
                    metric: float(row[metric].strip())
                    for metric in available_metrics
                }
            except (AttributeError, TypeError, ValueError) as error:
                raise ValueError(
                    f'line {line_number} has an invalid or missing required value'
                ) from error

            if offtarget_sequence != '':
                site = {
                    'seq': offtarget_sequence,
                    'pam': pam,
                    'disco_score': disco_score,
                    'coords': coords,
                    'mismatches': mismatches,
                }
                site.update(metric_values)
                offtargets.append(site)
    # Sort by disco score, then by mismatches
    offtargets.sort(
        key=lambda x: (-x['disco_score'], x['mismatches'], x['coords'])
    )
    return offtargets 


def visualizeOfftargets(infile: str, outfile: str, target_guide: str, pam: str, location: int,  drawseqs=True, title=None):

    # Get offtargets array from file
    offtargets  = parseSitesFile(infile)

    # Define top and left margins
    x_offset = 20
    y_offset = 50

    if location == 3:
        target_printseq = target_guide + pam
    elif location == 5:
        target_printseq = pam + target_guide
    else:
        raise ValueError("PAM location must be 5 or 3")

    available_metrics = list(offtargets.available_metrics)
    column_specs = [('disco_score', ('DISCO', 'Score'), 62)]
    column_specs.extend((metric, (metric,), 78 if metric == 'DPM' else 92)
                        for metric in available_metrics)
    column_specs.extend((
        ('mismatches', ('Mismatches',), 102),
        ('coords', ('Coordinates',), 230),
    ))
    column_positions = {}
    column_x = x_offset + box_size * len(target_printseq) + 8
    for key, unused_heading, width in column_specs:
        column_positions[key] = column_x
        column_x += width

    canvas_width = max(816, column_x + 20)
    canvas_height = max(1056, y_offset + box_size * (len(offtargets) + 3) + 20)
    dwg = svgwrite.Drawing(outfile, size=(canvas_width, canvas_height))

	# Header line for figure, use title if there is one, otherwise...bam file name?
    if title is None:
        dwg.add(dwg.text("DISCOVER-Seq Off Targets", insert=(x_offset, 30), style="font-size:20px; font-family:Arial"))
    else:
        # Define top and left margins
        dwg.add(dwg.text("DISCOVER-Seq Off Targets: "+title, insert=(x_offset, 30), style="font-size:20px; font-family:Arial"))

    x = 0
    if location == 5:
        x = box_size*(len(pam)-1) 
    else:
        x = x_offset+box_size*(len(target_guide))
    y = y_offset - 2
    dwg.add(dwg.text("PAM", insert=(x,y), style="font-size:14px; font-family:Courier"))
    # Draw reference sequence row
    for i, c in enumerate(target_printseq):
        y = y_offset
        if i < len(target_printseq):
            x = x_offset + i * box_size
        else:
            x = x_offset + i * box_size + box_size/5
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))

    y = y_offset + box_size + 5
    for key, heading_lines, unused_width in column_specs:
        for line_index, heading in enumerate(heading_lines):
            heading_y = y - 15 * (len(heading_lines) - line_index - 1)
            dwg.add(dwg.text(
                heading, insert=(column_positions[key], heading_y),
                style="font-size:15px; font-family:Arial"
            ))

    # Draw aligned sequence rows
    y_offset += 10  # leave gap after target sequence
    for j, seq in enumerate(offtargets):
        y = y_offset + j * box_size
        if location == 5:
            ot_printseq = seq['pam']+seq['seq']
        else:
            ot_printseq = seq['seq']+seq['pam']
        for i, c in enumerate(ot_printseq):
            if i < len(ot_printseq):
                x = x_offset + i * box_size
            else:
                x = x_offset + i * box_size + box_size/5
            if c == target_printseq[i] or target_printseq[i] == 'N':
                dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
            else:
                dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                if drawseqs:
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))

        y = y_offset + box_size * (j + 2) - 2
        sa = "font-size:14px; font-family:Arial"
        sc = "font-size:14px; font-family:courier"
        fill = 'red' if int(seq['mismatches']) == 0 else 'black'
        # If this is the target sequence, output it in red
        if int(seq['mismatches']) == 0:
            mismatch_value = "Target"
        else:
            mismatch_value = str(seq['mismatches'])
        values = {
            'disco_score': str(seq['disco_score']),
            'mismatches': mismatch_value,
            'coords': str(seq['coords']),
        }
        values.update({metric: f"{seq[metric]:.3f}"
                       for metric in available_metrics})
        for key, unused_heading, unused_width in column_specs:
            style = sc if key == 'coords' else sa
            dwg.add(dwg.text(
                values[key], insert=(column_positions[key], y),
                fill=fill, style=style
            ))
    dwg.save()

if __name__ == '__main__':
    args = parse_arguments()
    infile = args.infile
    outfile = args.outfile
    guide = args.guide
    pam = args.pam
    location = args.location
    title = args.title
    drawseqs = not args.notext
    visualizeOfftargets(infile, outfile, guide, pam, location, drawseqs, title=title)
