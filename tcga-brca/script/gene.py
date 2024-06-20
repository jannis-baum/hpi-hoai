class Gene:
    name: str
    id: str
    chromosome: str
    bp_start: int
    bp_end: int

    def __init__(self, name: str, id: str, chromosome: str, bp_start: int, bp_end: int):
        self.name = name
        self.id = id
        self.chromosome = chromosome
        self.bp_start = bp_start
        self.bp_end = bp_end
