from Bio import SeqIO
import csv

def find_sequence_positions(sequence, target):
    positions = []
    index = sequence.find(target)
    while index != -1:
        positions.append(index)
        index = sequence.find(target, index + 1)
    return positions

def find_distances(core_positions, ce_positions):
    distances = []
    for core_pos in core_positions:
        for ce_pos in ce_positions:
            distances.append(abs(core_pos - ce_pos))
    return distances

def process_sequence(record, core_sequence, ce_sequence):
    forward_core_positions = find_sequence_positions(record.seq, core_sequence)
    reverse_core_positions = find_sequence_positions(record.seq.reverse_complement(), core_sequence)

    forward_ce_positions = find_sequence_positions(record.seq, ce_sequence)
    reverse_ce_positions = find_sequence_positions(record.seq.reverse_complement(), ce_sequence)

    distances = find_distances(forward_core_positions, forward_ce_positions) + find_distances(reverse_core_positions, reverse_ce_positions)

    return forward_core_positions, reverse_core_positions, forward_ce_positions, reverse_ce_positions, distances

def main(input_file, core_sequence, ce_sequence, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['contig_name', 'core_sequence_orientation', 'core_sequence_start_position', 'core_sequence_end_position',
                      '37CE_sequence_orientation', '37CE_start_position', '37CE_end_position', 'distance']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for record in SeqIO.parse(input_file, 'fasta'):
            forward_core, reverse_core, forward_ce, reverse_ce, distances = process_sequence(record, core_sequence, ce_sequence)

            for core_pos in forward_core:
                for ce_pos in forward_ce:
                    writer.writerow({'contig_name': record.id,
                                     'core_sequence_orientation': 'forward',
                                     'core_sequence_start_position': core_pos,
                                     'core_sequence_end_position': core_pos + len(core_sequence),
                                     '37CE_sequence_orientation': 'forward',
                                     '37CE_start_position': ce_pos,
                                     '37CE_end_position': ce_pos + len(ce_sequence),
                                     'distance': abs(core_pos - ce_pos)})

            for core_pos in reverse_core:
                for ce_pos in reverse_ce:
                    writer.writerow({'contig_name': record.id,
                                     'core_sequence_orientation': 'reverse',
                                     'core_sequence_start_position': core_pos,
                                     'core_sequence_end_position': core_pos + len(core_sequence),
                                     '37CE_sequence_orientation': 'reverse',
                                     '37CE_start_position': ce_pos,
                                     '37CE_end_position': ce_pos + len(ce_sequence),
                                     'distance': abs(core_pos - ce_pos)})

if __name__ == "__main__":
    input_file = "/Users/input/cps_inter.fasta"
    core_sequence = "ttaccgtaaaaaagtga"
    ce_sequence = "ttgaaac"
    output_file = "output_results.csv"

    main(input_file, core_sequence, ce_sequence, output_file)
