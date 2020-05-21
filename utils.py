def chunk(iterable, size):
    """Yield successive n-sized chunks from iterable."""
    for i in range(0, len(iterable), size):
        yield iterable[i:i + size]


def group_ranges(values):
    """
    Chunk the input list into groups of consecutive elements.

    Each group is in form of [start_index, end_index)
    """
    units = []
    prev = values[0]

    for value in values:
        if value == prev + 1:
            units[-1].append(value)
            # None
        else:
            units.append([value])
        prev = value

    return [[u[0], u[-1] + 1] if len(u) > 1 else [u[0], u[0] + 1] for u in units]

def removeLn(file):
    import os, sys
    """
    remove last line file
    """
    readFile = open(file)
    lines = readFile.readlines()
    readFile.close()
    w = open(file,'w')
    w.writelines([item for item in lines[:-1]])
    w.close()
