import re
"""parser oggetti output clustal"""

from clustal import ClustalAlignment, ClustalSequence
import utils


class ClustalParser:
    def __init__(self, nseq=3):
        """
        Create a parser for a given number of sequences
        """
        self.nseq = nseq

    def parse(self, filename, reference=None):
        """
        Read the file in input assuming is a valid clustal omega file
        """
        with open(filename, 'r') as stream:
            sequences = self.parseLines(stream.readlines())

            return ClustalAlignment(
                reference=reference,
                sequences=sequences,
                length=max(int(s.bases) for s in sequences.values())
            )

    def parseLines(self, lines):
        """
        Parse each line in chunks of nseq + 2 items. Each block contains:
        - nseq sequence lines
        - one line for alignment outcome (skipped)
        - a separator line among blocks (skipped)
        """
        sequences = {}

        for block in utils.chunk(self._skipHeader(lines), self.nseq + 2):
            for line in block[0:self.nseq]:
                key, rawSequence, bases = self._parseLine(line)
                sequence = sequences.get(key, ClustalSequence(key)).extend(rawSequence, bases)

                if not key in sequences:
                    sequences[key] = sequence

        return sequences

    def _skipHeader(self, lines):
        """
        Skip the first n lines used by clustal to put comments
        """
        return lines[3:]

    def _parseLine(self, line):
        """
        Parse a single block line by stripping multiple consecutive space characters
        and replacing them with a single space, then split by the given space to obtain
        three parts.
        """
        return re.sub('\s+', ' ', line).strip(' ').split(' ')
