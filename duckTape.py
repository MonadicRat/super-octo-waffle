samples = open("sample.tsv")

samples.readline()
for sample in samples:
    fields = sample.split("\t")
    age = fields[1].replace("\"", "")
    age = age.split()[1][0:2]
    print(fields[0], age, sep = "\t")
