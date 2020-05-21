import utils


class ClustalAlignment:
    def __init__(self, reference, sequences=[], length=0):
        self.reference = reference
        self.sequences = sequences
        self.length = int(length)

    def __len__(self):
        return self.length

    def __getitem__(self, key):
        return self.peek_columns(key)

    def peek_reference(self, column, end=None):
        return self.peek_column(self.reference, column, end)

    def peek_others(self, column, end=None):
        return [{
            'sequence': sequence,
            'value': self.peek_column(sequence, column, end)
         } for sequence in self.sequences if sequence != self.reference]

    def peek_columns(self, column, end=None):
        return [self.peek_column(sequence, column, end) for sequence in self.sequences]

    def peek_column(self, key, column, end=None):
        end = end or column + 1
        return self.sequences[key][column:end]


class ClustalSequence:
    def __init__(self, id, sequence='', bases=None):
        self.id = id
        self.sequence = sequence
        self.bases = bases

    def extend(self, rawSequence, bases=None):
        self.sequence += rawSequence
        self.bases = bases

        return self

    def at(self, index):
        return self.sequence[index]

    def __len__(self):
        return self.bases if (self.bases != None) else len(self.sequence)

    def __repr__(self):
        return str(self)

    def __getitem__(self, key):
        return self.sequence[key]

    def __str__(self):
        return '<Sequence:{}, {}, {}>'.format(self.id, "self.sequence", self.bases)
