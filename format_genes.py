from string import printable

with open("genes_raw.txt", "r") as file:
    data = file.read().replace("\r\n", "\n")
    data = "".join(char for char in data if char in printable)
    data = data.replace(",", "\n")
    data = data.replace(" ", "")
    data = data.replace(";", "\n")
    f = open("genes.tsv", "w")
    f.write(data)
    f.close()
    Genes = set(line.strip() for line in open("genes.tsv"))
    contents = "\n".join(list(Genes))
    with open("genes.tsv", "w") as write:
        write.write("Gene\n" + contents)
