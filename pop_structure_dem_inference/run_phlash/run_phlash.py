import pandas as pd
import phlash
import os
import numpy as np
import argparse

#parse population argument
parser = argparse.ArgumentParser()
parser.add_argument('pop_arg', help='Population name')
args = parser.parse_args()

#load population (fam) file
fam = pd.read_csv("allpops.fam", delim_whitespace=True, header=None, names=["population", "sample"])
popmap = dict(zip(fam["sample"], fam["population"]))

#select Oahu samples....
oahu_samples = [sample for sample, pop in popmap.items() if args.pop_arg in pop]
print(oahu_samples)


#specify vcf file and regions
vcf_path = "repeat_filtered.vcf.gz"
scaffolds = ["scaffold_5", "scaffold_9", "scaffold_14"]

#load regions
contigs = []
for chrom in scaffolds:
    print(f"Loading {chrom}...")
    contigs.append(
        phlash.contig(
            vcf_path,
            samples=oahu_samples,
            region=f"{chrom}:10000000-130000000" 
        )
    )

#run phlash
results = phlash.fit(contigs)

times = np.array([dm.eta.t[1:] for dm in results])
#choose a grid of points at which to evaluate the size history functions
T = np.geomspace(times.min(), times.max(), 1000)
Nes = np.array([dm.eta(T, Ne=True) for dm in results])

#output
np.savetxt(f"Nes_phlash_{args.pop_arg}.csv", Nes, delimiter=",")
