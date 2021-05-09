
input_filename = r"data/cf1_S3_L001_R1_001.fastq"
output_filename = r"data/out.fastq"
num_reads = 20000

with open(input_filename, 'r') as in_file:
    with open(output_filename, "w") as out_file:
        for i in range(num_reads*4):
            line = in_file.readline()
            out_file.writelines("%s" % line)
            if i % 5000*4 == 0:
                print("{}/{} done.".format(i/4, num_reads))
