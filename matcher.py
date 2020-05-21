import utils


class AligmentDifferenceFinder:
    def analyze(self, alignment):
        self.unmatches = {}

        for sequence in alignment.sequences:
            if sequence != alignment.reference:
                self.unmatches[sequence] = []

                for i in range(alignment.length):
                    column = alignment.peek_column(sequence, i)

                    if column != 'N' and column != alignment.peek_reference(i):
                        self.unmatches[sequence].append(i)

                self.unmatches[sequence] = utils.group_ranges(self.unmatches[sequence])

        return self.unmatches
