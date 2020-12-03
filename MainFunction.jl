# remove the none coding framgents of the sequences
using DataFrames, CSV, Bio.Seq, ProgressMeter


cds = CSV.read("BarleyCDSvartest.txt", header=true, delim="\t")
# create a new dataframe where the number of mismatches and the count of stop codons are written - 8 new columns plus gene name
inN = DataFrame(Gene=String[], ISRStart=String[], GolfStart=String[], ISRStop=String[], GolfStop=String[], ISRStopcount=Int64[], GolfStopcount=Int64[], Mismatches=Int64[], DeletionISR=Int64[],  DeletionGolf=Int64[], Broken=String[], ProteinSize=Int64[], seqISR=String[], seqGolf=String[], seqCDS=String[])
# filter out the junk from the mitochondrien and chloroplast
cds = cds[occursin.("HORVU", cds.Gene),:]
sort!(cds, :Gene)

## test
#ay = ["HORVU1Hr1G000940","HORVU3Hr1G015510","HORVU7Hr1G000150","HORVU7Hr1G036390","HORVU7Hr1G073390","HORVU7Hr1G082540","HORVU1Hr1G000940","HORVU2Hr1G122620","HORVU3Hr1G015510","HORVU4Hr1G000020","HORVU4Hr1G072560","HORVU5Hr1G070050","HORVU6Hr1G009280","HORVU6Hr1G077620","HORVU6Hr1G084740","HORVU7Hr1G000150","HORVU7Hr1G036390","HORVU7Hr1G118560"]
#cds = cds[indexin(ay, cds.Gene)[.!(isnothing.(indexin(ay, cds.Gene)))],:] # check the copper genes
gen = "HORVU0Hr1G000430"
qwe2 = vcat(cds[cds.Gene.==gen,1], cds[cds.Gene.==gen,6])
loc = cds[cds.Gene.==gen,2][1]
dir = 0
###


cd("Protein_Comparisen") # set the working directory
# running a loop that extracts the Locus, gene and direction and transfering this to the Mega bash script
Threads.@threads for i in 1:size(cds,1)

        loc,gen = cds.Locus[i], cds.Gene[i]
        println(gen)
        cds.Direction[i] == "rev" ? dir = 1 : dir = 0
        run(`MegaX_gene_check_rev_cdscall.sh $loc $gen $dir`)
        # add the cds sequence to it
        qweg = CSV.read("golf.consensus$(gen)_2.fa", delim="\t", header=false)
        qweg = DataFrame(qweg)
        qwei = CSV.read("isr.consensus$(gen)_2.fa", delim="\t", header=false)
        qwei = DataFrame(qwei)

        if !(isempty(qweg)) && !(isempty(qwei)) # when the sequence of one genotype cannot be called
            qweg[1,1]=">Golf"
            qwei[1,1]=">ISR42-8"
            qwe2 = vcat(cds.SpliceVaraint[i],cds.Seq[i])
            qwe = DataFrame(vcat(Matrix(qweg), Matrix(qwei),qwe2))
            CSV.write("MegaX_$gen.fas", qwe, delim="\t",  writeheader=false)

            #get the dcs aligned to the genomic dna
            run(`megacc -a uscle_align_nucleotide.mao -d MegaX_$gen.fas -f Fasta -o $gen`)
            # make a single line file out of this
            run(`bash awk_helper.sh $gen`)
            # remove the first blank line
            run(`sed -i 1d $(gen)_2.fasta`)
            # push the data into two columns
            #run(`xargs -n2 < /home/michael/Schreibtisch/TEst/test2.fasta > /home/michael/Schreibtisch/TEst/test4.fasta`)



            ###################################### Julia code ######################################
            seq = CSV.read("$(gen)_2.fasta", header=false, delim="\t")
            seq = seq[[2,4,6],1]

            seq = hcat(["Golf", "ISR42-8", "CDS"], seq)

            # get the positions where it is matched
            x = findall(r"A|C|T|G|N", seq[3,2])
            y = Int64[]
            for i in 1:length(x); append!(y,parse(Int64,split.(string.(x),":")[i][1])); end

            # extract these from the seq of all 4 genotypes
            for i in 1:3; seq[i,2] = seq[i,2][y] ; end

            # get the amino acid information

            cdS = translate(DNASequence(seq[3,2]))
                # first check - is the seq propper aligned? -> starting with ATG, ending with TAA|TGA|TAG
                # check if the stop codon is there


                if occursin("-",seq[1,2]) == false
                                                    golf = translate(DNASequence(seq[1,2]))
                                                    delG = 0
                else
                        #if startswith.(seq[2,2], "ATG") && endswith.(seq[2,2], r"TAA|TGA|TAG") # Golf
                                                                                                        # remove the gaps for Golf
                                                                                                        # get the number of deletions
                                                                                                        delG =  length(findall(r"-", seq[1,2]))
                                                                                                        x = findall(r"A|C|T|G|N", seq[1,2])
                                                                                                        y = Int64[]
                                                                                                        for i in 1:length(x); append!(y,parse(Int64,split.(string.(x),":")[i][1])); end
                                                                                                        seq[1,2] = seq[1,2][y]

                                                                                                        # only if the sequence is a multiple of 3, the translation will be done - other wise it will cause an error
                                                                                                        # if the seq has a deletion or insertion of 1,2,4,5,7,etc the DNA can not be translated - to avoid this, the sequence should be completed as with what might be missing
                                                                                                        if mod(length(seq[1,2]),3) == 0
                                                                                                                                            golf = translate(DNASequence(seq[1,2]))
                                                                                                        elseif mod(length(seq[1,2]),3) == 2 # one base missing
                                                                                                                                            golf = seq[1,2] * SubString(seq[3,2],length(seq[3,2]))
                                                                                                                                            golf = translate(DNASequence(golf))
                                                                                                        elseif mod(length(seq[1,2]),3) == 1 # one base missing
                                                                                                                                            golf = seq[1,2] * SubString(seq[3,2],length(seq[3,2])-1)
                                                                                                                                            golf = translate(DNASequence(golf))
                                                                                                        else
                                                                                                                                            golf = 1
                                                                                                        end
                        #else
                        #    golf = 1

                        #end
                end

                if occursin("-",seq[2,2]) == false
                                                    isr = translate(DNASequence(seq[2,2]))
                                                    delI = 0
                else

                        #if startswith.(seq[3,2], "ATG") && endswith.(seq[3,2], r"TAA|TGA|TAG") # ISR 42-8
                                                                                                # get the number of deletions
                                                                                                delI =  length(findall(r"-", seq[2,2]))
                                                                                                x = findall(r"A|C|T|G|N", seq[2,2])
                                                                                                y = Int64[]
                                                                                                for i in 1:length(x); append!(y,parse(Int64,split.(string.(x),":")[i][1])); end
                                                                                                seq[2,2] = seq[2,2][y]

                                                                                                if mod(length(seq[2,2]),3) == 0
                                                                                                                                    isr = translate(DNASequence(seq[2,2]))
                                                                                                elseif mod(length(seq[2,2]),3) == 2 # one base missing
                                                                                                                                    isr = seq[2,2] * SubString(seq[3,2],length(seq[3,2]))
                                                                                                                                    isr = translate(DNASequence(isr))
                                                                                                elseif mod(length(seq[2,2]),3) == 1 # one base missing
                                                                                                                                    isr = seq[2,2] * SubString(seq[3,2],length(seq[3,2])-1)
                                                                                                                                    isr = translate(DNASequence(isr))
                                                                                                else
                                                                                                                                    isr = 1
                                                                                                end
                      #else
                        #  isr = 1
                        #end

                end
            # check for stop codons in the proteine
            csg = count("*",String(golf))
            csi = count("*",String(isr))

            # check if Golf and / or ISR have been translated
            if isr != 1 && golf != 1 # both okay
                                    # how many AS are differnt?
                                    qay = Int64[]
                                    # check which is the shortest and take this value
                                    vb = minimum([length(cdS),length(golf), length(isr)])
                                    for i in 1:vb ; golf[i] == isr[i] ? append!(qay,1) : append!(qay,0);end
                                    info = sum(qay.==0)

            else # if one or both are broken
                                    info = 999
            end

            push!(inN,[gen, string(startswith.(seq[2,2], "ATG")), string(startswith.(seq[1,2], "ATG")), string(endswith.(seq[2,2], r"TAA|TGA|TAG")), string(endswith.(seq[1,2], r"TAA|TGA|TAG")),
            csi, csg, info, delI, delG, "no", length(isr), String(isr), String(golf), String(cdS)])
            # write to table, so that if an error occurse at the end of the loop not all data is lost
            run(` rm golf.consensus$(gen)_2.fa isr.consensus$(gen)_2.fa $(gen)_summary.txt $gen.fasta  MegaX_$gen.fas`)

        else
            isempty(qwei) ? broke = "ISR42-8" : broke = "Golf"
            push!(inN,[gen,"","","","",0,0,0,0,0, broke, 0, "", "" ,""])
        end

            println(last(inN[!,1:12],1))
            # CSV.write("/Protein_comparisen.txt", inN; delim="\t", append=true)

end

CSV.write("Protein_comparisenALL.txt", inN; delim="\t")
