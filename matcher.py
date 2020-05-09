class AligmentDifferenceFinder:
    def __init__(self):
        self.unmatches = []
        pass

    def analyze(self, alignment):
        self.unmatches = []

        for i in range(alignment.length):
            column = alignment[i]

            j = 0
            while j < len(column) and column[j] == column[(j+1) % len(column)]:
                j += 1

            if j != len(column):
                self.unmatches.append(i)

        # return utils.group_ranges(self.unmatches)
        return self.unmatches
