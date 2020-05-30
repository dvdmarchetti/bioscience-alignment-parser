import utils


class AligmentDifferenceFinder:
    def analyze(self, alignment):
        self.unmatches = {}

        for sequence in alignment.sequences:
            if sequence != alignment.reference:
                self.unmatches[sequence] = []
                #iterate alignment
                for i in range(alignment.length):
                    column = alignment.peek_column(sequence, i)

                    ref_char = alignment.peek_reference(i)
                    if (
                        (column == 'N') or
                        (column == 'R' and (ref_char == 'A' or ref_char == 'G')) or
                        (column == 'Y' and (ref_char == 'C' or ref_char == 'T')) or
                        (column == 'S' and (ref_char == 'G' or ref_char == 'C')) or
                        (column == 'W' and (ref_char == 'A' or ref_char == 'T')) or
                        (column == 'K' and (ref_char == 'G' or ref_char == 'T')) or
                        (column == 'M' and (ref_char == 'A' or ref_char == 'C')) or
                        (column == 'B' and (ref_char == 'C' or ref_char == 'G' or ref_char == 'T')) or
                        (column == 'D' and (ref_char == 'A' or ref_char == 'G' or ref_char == 'T')) or
                        (column == 'H' and (ref_char == 'A' or ref_char == 'C' or ref_char == 'T')) or
                        (column == 'V' and (ref_char == 'A' or ref_char == 'C' or ref_char == 'G'))
                    ):
                        continue

                    if column != ref_char:
                        self.unmatches[sequence].append(i)

                self.unmatches[sequence] = utils.group_ranges(self.unmatches[sequence])

        return self.unmatches
