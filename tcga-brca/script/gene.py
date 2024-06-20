class Gene:
    name: str
    chromosome: str
    bp_start: int
    bp_end: int

    def __init__(self, name, chromosome, bp_start, bp_end):
        self.name = name
        self.chromosome = chromosome
        self.bp_start = bp_start
        self.bp_end = bp_end
