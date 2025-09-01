from random import choice, random, randint
import pandas 
from tqdm import tqdm
alphabet = set("ACGT")

class Genome:
    def __init__(self, sequence = 1000, parent = None):
        
        if type(sequence) is int:
            _alphabet = list(alphabet)
            sequence = "".join(choice(_alphabet) for _ in range(sequence))

        self.sequence = sequence
        self.parent = parent

    def __repr__(self):
        return f"Genome(sequence=[{alphabet}]*{len(self.sequence)})"

    def mutate(self, rate = 0.001, inplace = False):

        new_sequence = "".join([i if random() > rate else choice(list(alphabet - {i})) for i in self.sequence])

        if inplace:
            self.sequence = new_sequence
        else:
            return Genome(sequence=new_sequence, parent=self)

    def __eq__(self, other):
        if not isinstance(other, Genome):
            return NotImplemented
        return all([a == b for a,b in zip(self.sequence, other.sequence)]) and self.parent == other.parent

    def homologous_recomb(self, other, length):
        if not isinstance(other, Genome):
            return NotImplemented
        start = randint(0, len(self.sequence) - length)
        end = start + length
        new_sequence = self.sequence[:start] + other.sequence[start:end] + self.sequence[end:]
        assert len(self.sequence) == len(other.sequence), "Recombination must not change the length of the genome"
        self.sequence = new_sequence


    def print_ali(self, other):
        if not isinstance(other, Genome):
            return NotImplemented
        assert len(self.sequence) == len(other.sequence), "Genomes must be of the same length to print alignment"
        sames = "".join(["." if a == b else " " for a, b in zip(self.sequence, other.sequence)])
        print(f"{self.sequence}\n{sames}\n{other.sequence}")

    def ani(self, other):
        if not isinstance(other, Genome):
            return NotImplemented
        matches = sum(a == b for a, b in zip(self.sequence, other.sequence))
        return (matches / len(self.sequence)) * 100

    @classmethod
    def pure_drift(cls, pop_size=100, generations=100, mutation_rate=0.001):

        luca = Genome()
        pop_size = pop_size
        pop = [luca.mutate() for i in range(pop_size)]

        for gen in tqdm(range(generations)):
            _ = [p.mutate(mutation_rate, inplace=True) for p in pop]

        all_anis = [ p.ani(q) for i,p in enumerate(pop) for j,q in enumerate(pop) if i > j]

        data = pandas.DataFrame.from_records({ 'anis' : all_anis, 'model' : 'pure_drift', 'pop_size' : pop_size, 'mutation_rate' : mutation_rate, 'recomb_rate' : [0]})
