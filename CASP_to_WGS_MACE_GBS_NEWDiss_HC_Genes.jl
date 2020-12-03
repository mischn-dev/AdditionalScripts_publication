using DelimitedFiles, DataFrames, Statistics, CSV, ProgressMeter, HypothesisTests, GLM, Crayons, Gadfly, Compose, CuArrays, RCall
using Cairo

# rato  = CASP ISR42-8 AlleleFrequency column
casp = CSV.read("/media/michael/Uni/Allelfreq/Pool_vergleich_GBS_MACE_WGS/Kaspar_auswertungGeneAnnotation3.csv")
ca = casp[:,[4,6,7,8,19,20,21,22]]
deleterows!(ca,(1,5))
ca[:,:CP] = string.(ca[:,:Chr],"_", ca[:,:Pos])

# relate casp to marker and to cintig positions

freq = readDict("Allelefrequency_EachGeneration_Marker", "std")

function CP(x)
              x[:,:CP] = string.(x[:,:Chr],"_", x[:,:Pos])
              y = join(ca, x, on=:CP, makeunique=true)
              println(cor(y[:,:ratio], y[:,:Allelfreq_Isr]))
              return y
 end  # every snp gets a unique identifer to match with the casp marker

function HP2(x)
                y= join(ca, x, on=:Gene_ID)
                 println(cor(y[:,:ratio], y[:,:Allelfreq_Isr]))
                 return y
    end # join the gene haplotypes

function HP(x)
            y = join(x, ca, on=(:Gene_ID, :MARKER), makeunique=true)
            println(cor(y[:,:ratio], y[:,:Allelfreq_Isr]))
            return y
      end # join the marker

function HCP(x)
            y = join(x, ca, on=(:Gene_ID, :Contig), makeunique=true)
            println(cor(y[:,:ratio], y[:,:Allelfreq_Isr]))
            return y
      end # join the contigs

function HPh(x)
           y = join(x, ca, on=(:Gene_ID, :Gene_ID_HC), makeunique=true)
           println(cor(y[:,:ratio], y[:,:Allelfreq_Isr]))
         return y
   end # join the high confidence genes


contig(freq, "genoytpes_wgs.txt")
# single snp comparisen
eps = CP(freq[:F22P1KEP]) # 0.9328779263442215 for 10 snps
CP(freq[:F23P1K1]) # 0.7519415112708296
CP(freq[:F23P1Ö1]) # 0.4518162283425448
CP(freq[:F22P1Ö1]) # 0.24538434885207874
CP(freq[:F3P1K1]) # -0.3555570972262029
CP(freq[:F16P1K1]) # 0.7386285068405248

# genetic map haplotypes
epg = HP(freq[:F22P1KEP_Haplotypes]) # 17 marker , 0.9579173953800624
HP(freq[:F23P1K1_Haplotypes]) # 0.9358692077740327
HP(freq[:F23P1Ö1_Haplotypes]) # 0.6229479186525863
HP(freq[:F22P1Ö1_Haplotypes]) # 0.7651550788391824
HP(freq[:F3P1K1_Haplotypes]) # 0.3597842483479961
HP(freq[:F16P1K1_Haplotypes]) # 0.9327137172855968

freq2 = readDict("Allelefrequency_EachGeneration_Genes", "std")

# physsical map haplotypes
ep = HP2(freq2[:F22P1KEP_Haplotypes]) # 0.9009249872725319
HP2(freq2[:F23P1K1_Haplotypes]) # 0.8989024415186938
HP2(freq2[:F23P1Ö1_Haplotypes]) # 0.5462944745261715
HP2(freq2[:F22P1Ö1_Haplotypes]) # 0.8176122864697924
HP2(freq2[:F3P1K1_Haplotypes]) # 0.39757522105094
HP2(freq2[:F16P1K1_Haplotypes]) # 0.9261194498703404

# contigs
epc = HCP(freq[:F22P1KEP_Contigs]) # 19 contigs, 0.9459697932748631
HCP(freq[:F22P1K2_Contigs])  # 19 contigs, 0.9263071927796536


# read coverage per sample average
median(eps.Readcount) # 7
median(epg.Readcount) # 7317
median(ep.Readcount) # 197
median(epc.Readcount) # 75461

# root mean square error RMSE
sqrt(sum(abs.(eps.Allelfreq_Isr .- eps.ratio).^2)/size(eps,1)) # single SNP 0.029303135072984238
sqrt(sum(abs.(epg.Allelfreq_Isr .- epg.ratio).^2)/size(epg,1)) # marker     0.007655437656637487
sqrt(sum(abs.(ep.Allelfreq_Isr .- ep.ratio).^2)/size(ep,1)) # gen      0.011032785456895584
sqrt(sum(abs.(epc.Allelfreq_Isr .- epc.ratio).^2)/size(epc,1)) # contig       0.008131288825251215


# zero inflated negative binomial model testing the casp to pool data
wgs = Dict(:snp => eps, :marker => epg, :gene => ep, :contig => epc)
for i in collect(keys(wgs))
            println(i)
            m = vcat(repeat([1], size(wgs[i],1)), repeat([2], size(wgs[i],1)))
            n = vcat(wgs[i][:,:Allelfreq_Isr], wgs[i][:,:ratio])
            n = Int64.(floor.(n .* 100))
            mn = DataFrame(hcat(m,n))

            if size(mn,1) > 2
            @rput mn
            R"""
            library(MASS)
            library(pscl)
            mn[,1] = as.factor(mn[,1])
            m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
            print(m1)
            """
         end
end


# filter - mimimum coverage for the haplotypes = 100 Reads
epg = epg[epg[!,:Readcount].>100,:]
ep = ep[ep[!,:Readcount].>100,:]


# plot ep # physical map haplotyes
p1 = plot(ep, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, ymin=:MinFreq, ymax=:MaxFreq, Geom.point, Geom.errorbar),
   layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
   layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
    Guide.xlabel(""), Guide.ylabel(""), Guide.colorkey(title="", pos=[0.7w,0.001h]),
    #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
    Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "B"))),Coord.cartesian(ymin=0.0,ymax=1.0),
    Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
    key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))
# plot ep single snps
p2 =  plot(eps, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, Geom.point),
   layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
   layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
    Guide.xlabel(""), Guide.ylabel("Pool Allele Frequency"), Guide.colorkey(title="", pos=[0.7w,0.001h]) ,
    #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
    Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "A"))),Coord.cartesian(ymin=0.0,ymax=1.0),
    Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
    key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))

# plot ep # genetic map
p3 =   plot(epg, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, ymin=:MinFreq, ymax=:MaxFreq, Geom.point, Geom.errorbar),
   layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
   layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
    Guide.xlabel("KASPar Allele Frequency"), Guide.ylabel("Pool Allele Frequency"), Guide.colorkey(title="", pos=[0.7w,0.001h]) ,
    #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
    Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "C"))),Coord.cartesian(ymin=0.0,ymax=1.0),
    Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
    key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))

p4 =   plot(epc, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, ymin=:MinFreq, ymax=:MaxFreq, Geom.point, Geom.errorbar),
       layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
       layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
        Guide.xlabel("KASPar Allele Frequency"), Guide.ylabel(""), Guide.colorkey(title="Coverage", pos=[0.7w,0.001h]) ,
        #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
        Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "D"))),Coord.cartesian(ymin=0.0,ymax=1.0),
        Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
        key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))

myplot = hstack(vstack(p2, p3), vstack(p1,p4))

draw(PNG("Casp_to_WGS_HC.png", 16inch, 9inch, dpi = 300),myplot)
####

## GBS correlation to CASP

wgs = "LGC_samtoolsGBS_lgc_snpcalling.vcf"
genolist = "LGC_genotypes.txt"
info = "160517_Hv_IBSC_PGSB_r1_HC_functional_annotation_sort.txt"
inf = "Marker_reference_all.txt"

feg = Allelfreqcalling(wgs, genolist)
genloc(inf, feg, genolist, true)
haplotyping(feg, genolist, "Marker")
contig(feg, genolist)

mg2 = CP(feg[:f22epkd]) # 1 snp
CP(feg[:f22iikd]) # 1 snp
CP(feg[:f22iiikd]) # 1 snp

mg3 = HP(feg[:f22epkd_Haplotypes]) # 19 marker - 0.834
HP(feg[:f22iikd_Haplotypes]) # 19 m - 0.759
HP(feg[:f22iiikd_Haplotypes]) # 19 m - 0.876

mg4 = HCP(feg[:f22epkd_Contigs]) # 19m - 0.944
HCP(feg[:f22iikd_Contigs]) # 19m - 0.925
HCP(feg[:f22iiikd_Contigs]) # 19m - 0.942

fep = Allelfreqcalling(wgs, genolist)
genloc(info, fep, genolist)
haplotyping(fep, genolist)

mg1 = HPh(fep[:f22epkd_Haplotypes]) # 9 gene - -0.06
HPh(fep[:f22iikd_Haplotypes]) # 9 gene - 0.773
HPh(fep[:f22iiikd_Haplotypes]) # 9 genes - 0.769

# read coverage per sample average
median(mg1.Readcount) # 56
median(mg2.Readcount) # 41
median(mg3.Readcount) # 575
median(mg4.Readcount) # 10122

# root mean square error RMSE
sqrt(sum(abs.(mg2.Allelfreq_Isr .- mg2.ratio).^2)/19) # single SNP  0.04172764396634686
sqrt(sum(abs.(mg3.Allelfreq_Isr .- mg3.ratio).^2)/17) # marker      0.0368150766259218
sqrt(sum(abs.(mg4.Allelfreq_Isr .- mg4.ratio).^2)/19) # contig      0.029643121128601957
sqrt(sum(abs.(mg1.Allelfreq_Isr .- mg1.ratio).^2)/19) # gene        0.037647295465793804

# zero inflated negative binomial model testing the casp to pool data
gbs = Dict(:snp => mg2, :marker => mg3, :gene => mg1, :contig => mg4)
for i in collect(keys(gbs))
            println(i)
            m = vcat(repeat([1], size(gbs[i],1)), repeat([2], size(gbs[i],1)))
            n = vcat(gbs[i][:,:Allelfreq_Isr], gbs[i][:,:ratio])
            n = Int64.(floor.(n .* 100))
            mn = DataFrame(hcat(m,n))

            if size(mn,1) > 2
            @rput mn
            R"""
            library(MASS)
            library(pscl)
            mn[,1] = as.factor(mn[,1])
            m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
            print(m1)
            """
         end
end

# plot ep gene haplotypes
g1 = plot(mg1, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, ymin=:MinFreq, ymax=:MaxFreq, Geom.point, Geom.errorbar),
   layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
   layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
    Guide.xlabel(""), Guide.ylabel(""),Guide.colorkey(title="", pos=[0.7w,0.001h]) ,
    #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
    Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "B"))),Coord.cartesian(ymin=0.0,ymax=1.0),
    Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
    key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))

# plot ep single snps
g2 =  plot(mg2, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, Geom.point),
   #layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
   layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
    Guide.xlabel(""), Guide.ylabel("Pool Allele Frequency"), Guide.colorkey(title="", pos=[0.7w,0.001h]) ,
    #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
    Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "A"))),Coord.cartesian(ymin=0.0,ymax=1.0),
    Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
    key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))

# plot ep # genetic map
g3 =   plot(mg3, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, ymin=:MinFreq, ymax=:MaxFreq, Geom.point, Geom.errorbar),
   layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
   layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
    Guide.xlabel("KASPar Allele Frequency"), Guide.ylabel("Pool Allele Frequency"), Guide.colorkey(title="", pos=[0.7w,0.001h]) ,
    #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
    Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "C"))),Coord.cartesian(ymin=0.0,ymax=1.0),
    Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
    key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))

# plot ep on contig level
g4 =   plot(mg4, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, ymin=:MinFreq, ymax=:MaxFreq, Geom.point, Geom.errorbar),
   layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
   layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
    Guide.xlabel("KASPar Allele Frequency"), Guide.ylabel(""), Guide.colorkey(title="Coverage", pos=[0.7w,0.001h]) ,
    #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
    Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "D"))),Coord.cartesian(ymin=0.0,ymax=1.0),
    Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
    key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))


GBSplot = hstack(vstack(g2, g3), vstack(g1, g4))
draw(PNG("DFGCasp_to_GBS_HC.png", 16inch, 9inch, dpi = 300), GBSplot)



#### MACE correlation to CASP - using genloc and haplotyping from version 4.7 - 6.8 is not assigning the haplotypes correclty
genomace = "MACE_DupRmf8f21f22.vcf"
genomg = "MACE_dup_genotypes.txt"

fepm = Allelfreqcalling(genomace, genomg)
genloc(info, fepm, genomg)
haplotyping(fepm, genomg)
contig(fepm, genomg)

m1 = HPh(fepm[:f22epkd_Haplotypes]) # 0.8837476002109095
HPh(fepm[:f22iikd_Haplotypes]) # 0.8647987694836098
HPh(fepm[:f22iiikd_Haplotypes]) # 0.8761137123110734

m2 = CP(fepm[:f22epkd]) # 0.7885271086614456 - 19 snps
CP(fepm[:f22iikd]) # 0.766971647051518
CP(fepm[:f22iiikd]) # 0.7295967808983826

m4 = HCP(fepm[:f22epkd_Contigs]) # 19m - 0.8854526323568372
HCP(fepm[:f22iikd_Contigs]) # 19m - 0.8482250438137181
HCP(fepm[:f22iiikd_Contigs]) # 19m - 0.8960142189162377

fegm = Allelfreqcalling(genomace, genomg)
genloc(inf, fegm, genomg, true)
haplotyping(fegm, genomg, "Marker")

m3 = HP(fegm[:f22epkd_Haplotypes]) # 0.9332791088253864
HP(fegm[:f22iikd_Haplotypes]) # 0.9244808237144871
HP(fegm[:f22iiikd_Haplotypes]) # 0.9709821114024797

# read coverage per sample average
median(m1.Readcount) # 66
median(m2.Readcount) # 51
median(m3.Readcount) # 158
median(m4.Readcount) # 917

# root mean square error RMSE
sqrt(sum(abs.(m2.Allelfreq_Isr .- m2.ratio).^2)/19) # single SNP  0.04172764396634686
sqrt(sum(abs.(m3.Allelfreq_Isr .- m3.ratio).^2)/17) # marker      0.0368150766259218
sqrt(sum(abs.(m4.Allelfreq_Isr .- m4.ratio).^2)/19) # contig      0.029643121128601957
sqrt(sum(abs.(m1.Allelfreq_Isr .- m1.ratio).^2)/19) # gene        0.037647295465793804

# zero inflated negative binomial model testing the casp to pool data
mace = Dict(:snp => m2, :marker => m3, :gene => m1, :contig => m4)
for i in collect(keys(mace))
            println(i)
            m = vcat(repeat([1], size(mace[i],1)), repeat([2], size(mace[i],1)))
            n = vcat(mace[i][:,:Allelfreq_Isr], mace[i][:,:ratio])
            n = Int64.(floor.(n .* 100))
            mn = DataFrame(hcat(m,n))

            if size(mn,1) > 2
            @rput mn
            R"""
            library(MASS)
            library(pscl)
            mn[,1] = as.factor(mn[,1])
            m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
            print(m1)
            """
         end
end


# plot
# plot ep # physical map haplotyes
pm1 = plot(m1, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, ymin=:MinFreq, ymax=:MaxFreq, Geom.point, Geom.errorbar),
   layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
   layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
    Guide.xlabel(""), Guide.ylabel(""), Guide.colorkey(title="", pos=[0.7w,0.001h]) ,
    #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
    Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "B"))),Coord.cartesian(ymin=0.0,ymax=1.0),
    Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
    key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))

# plot ep single snps
pm2 =  plot(m2, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, Geom.point),
   layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
   layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
    Guide.xlabel(""), Guide.ylabel("Pool Allele Frequency"), Guide.colorkey(title="", pos=[0.7w,0.001h]) ,
    #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
    Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "A"))),Coord.cartesian(ymin=0.0,ymax=1.0),
    Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
    key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))
# plot ep # genetic map
pm3 =   plot(m3, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, ymin=:MinFreq, ymax=:MaxFreq, Geom.point, Geom.errorbar),
   layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
   layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
    Guide.xlabel("KASPar Allele Frequency"), Guide.ylabel("Pool Allele Frequency"), Guide.colorkey(title="", pos=[0.7w,0.001h]) ,
    #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
    Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "C"))),Coord.cartesian(ymin=0.0,ymax=1.0),
    Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
    key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))

pm4 =   plot(m4, layer(x=:ratio, y=:Allelfreq_Isr, color=:Readcount, ymin=:MinFreq, ymax=:MaxFreq, Geom.point, Geom.errorbar),
       layer(x=:ratio, y=:Allelfreq_Isr, Geom.smooth, Theme(default_color=colorant"firebrick")),
       layer(x=collect(Float64, 0:0.1:0.4), y= collect(Float64, 0:0.1:0.4), Geom.abline(color="gray", style=:dash)),
        Guide.xlabel("KASPar Allele Frequency"), Guide.ylabel(""), Guide.colorkey(title="Coverage", pos=[0.7w,0.001h]) ,
        #Scale.color_discrete_manual("sienna4","sienna3","sienna2","sienna1", "tan1", "peachpuff2", order=[5,4,6,3,2,1]),
        Guide.annotation(compose(context(), Compose.text(0.01, 0.9, "D"))),Coord.cartesian(ymin=0.0,ymax=1.0),
        Theme(background_color=colorant"white", default_color="black", major_label_color=colorant"black", minor_label_color=colorant"black", major_label_font_size=16pt,minor_label_font_size=14pt,
        key_title_font_size=16pt, key_title_color=colorant"black", key_label_font_size=14pt, key_label_color=colorant"black"))


mp2 = hstack(vstack(pm2,pm3), vstack(pm1,pm4))
draw(PNG("MACEtoCASP_HC.png", 16inch, 9inch, dpi=300), mp2)

pp = hstack(myplot, mp2, GBSplot)
draw(PNG("CASP_vs_MACE_vs_WGS_vs_GBS_letter.png", 20inch, 11inch,  dpi=500), pp)


############# compare mace vs WGS for the EP dataset #############################################################

# for the casp marker
wgs = HP2(freq2[:F22P1KEP_Haplotypes])
mace = HP2(fepm[:f22epkd_Haplotypes])
wm = join(wgs, mace, on=:Gene_ID , makeunique=true)
cor(wm[:,:Allelfreq_Isr], wm[:,:Allelfreq_Isr_1]) # 0.765 for 19 Marker

## for the whole genome

# snp level

function CP(x,y)
              x[:,:CP] = string.(x[:,:Chr],"_", x[:,:Pos])
              y[:,:CP] = string.(y[:,:Chr],"_", y[:,:Pos])
              z = join(x, y, on=:CP, makeunique=true)
              z = z[.&(z[:,:Alleles_1].!= "./.", z[:,:Alleles].!= "./."),:]
              println(cor(z[:,:Allelfreq_Isr], z[:,:Allelfreq_Isr_1]))
              return z
 end

z = CP(freq2[:F22P1KEP], fep[:f22epkd])

z[:,:Allelfreq_Isr_1] = Float64.(z[:,:Allelfreq_Isr_1])
z = z[.!(isnan.(z[:,:Allelfreq_Isr_1])),:]
cor(z[:,:Allelfreq_Isr], z[:,:Allelfreq_Isr_1]) # 0.45183542259366777
@rput z
R"""
plot(z$Allelfreq_Isr, z$Allelfreq_Isr_1)
"""

# get the information of the statistics how to interpret the variation
m = vcat(repeat([1], size(z,1)), repeat([2], size(z,1)))
n = vcat(z[:,:Allelfreq_Isr], z[:,:Allelfreq_Isr_1])
n = Int64.(floor.(n .* 100))
mn = DataFrame(hcat(m,n))

@rput mn
R"""
library(MASS)
library(pscl)
mn[,1] = as.factor(mn[,1])
m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]
"""
# negbino 8.430902e-09
# zeroinfl 4.324915e-20


#####################################################
## haplotype level
######################################################
# genetic map MACE to WGS

gen = join(fegm[:f22epkd_Haplotypes], freq[:F22P1KEP_Haplotypes], on=:Gene_ID, makeunique=true)
cor(gen[:,:Allelfreq_Isr], gen[:,:Allelfreq_Isr_1]) # 0.6523240087293424 for 3067 Haplotypes

# get the information of the statistics how to interpret the variation
m = vcat(repeat([1], size(gen,1)), repeat([2], size(gen,1)))
n = vcat(gen[:,:Allelfreq_Isr], gen[:,:Allelfreq_Isr_1])
n = Int64.(floor.(n .* 100))
mn = DataFrame(hcat(m,n))

@rput mn
R"""
library(MASS)
library(pscl)
mn[,1] = as.factor(mn[,1])
m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]
"""
# 3.977102e-40
# zeroinfl   0.36

@rput gen
R"""
plot(gen$Allelfreq_Isr, gen$Allelfreq_Isr_1)
"""

##################################################
# physikalische karte haplotypen MACE to WGS
phy = join(fepm[:f22epkd_Haplotypes], freq2[:F22P1KEP_Haplotypes], on=:Gene_ID, makeunique=true)
cor(phy[:,:Allelfreq_Isr], phy[:,:Allelfreq_Isr_1]) # 0.5767062992568456 of 5377 Haplotypes

# get the information of the statistics how to interpret the variation
m = vcat(repeat([1], size(phy,1)), repeat([2], size(phy,1)))
n = vcat(phy[:,:Allelfreq_Isr], phy[:,:Allelfreq_Isr_1])
n = Int64.(floor.(n .* 100))
mn = DataFrame(hcat(m,n))

@rput mn
R"""
library(MASS)
library(pscl)
mn[,1] = as.factor(mn[,1])
m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]
"""
# negbino 8.830953e-26
# zeroinfl 5.274772e-01

@rput phy
R"""
plot(phy$Allelfreq_Isr, phy$Allelfreq_Isr_1)
"""

################################################## WGS to GBS ###
# gene level

phy = join(fep[:f22epkd_Haplotypes], freq2[:F22P1KEP_Haplotypes], on=:Gene_ID, makeunique=true)
cor(phy[:,:Allelfreq_Isr], phy[:,:Allelfreq_Isr_1]) # 0.010600093122863667 of 6107 Haplotypes

# get the information of the statistics how to interpret the variation
m = vcat(repeat([1], size(phy,1)), repeat([2], size(phy,1)))
n = vcat(phy[:,:Allelfreq_Isr], phy[:,:Allelfreq_Isr_1])
n = Int64.(floor.(n .* 100))
mn = DataFrame(hcat(m,n))

@rput mn
R"""
library(MASS)
library(pscl)
mn[,1] = as.factor(mn[,1])
m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]
"""
# negbino <2e-16
# zeroinfl <2e-16

# marker level
gen = join(feg[:f22epkd_Haplotypes], freq[:F22P1KEP_Haplotypes], on=:Gene_ID, makeunique=true)
cor(gen[:,:Allelfreq_Isr], gen[:,:Allelfreq_Isr_1]) # 0.6938072142057771 of 3686 Haplotypes

# get the information of the statistics how to interpret the variation
m = vcat(repeat([1], size(gen,1)), repeat([2], size(gen,1)))
n = vcat(gen[:,:Allelfreq_Isr], gen[:,:Allelfreq_Isr_1])
n = Int64.(floor.(n .* 100))
mn = DataFrame(hcat(m,n))

@rput mn
R"""
library(MASS)
library(pscl)
mn[,1] = as.factor(mn[,1])
m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]
"""
# negbino <2e-16
# zeroinfl 0.502



#### casp vs WGS comaprisen with stats test ###############################################################

### EP - gene haplotypes

# get the information of the statistics how to interpret the variation
m = vcat(repeat([1], size(ep,1)), repeat([2], size(ep,1)))
n = vcat(ep[:,:ratio], ep[:,:Allelfreq_Isr])
n = Int64.(floor.(n .* 100))
mn = DataFrame(hcat(m,n))

@rput mn
R"""
library(MASS)
library(pscl)
mn[,1] = as.factor(mn[,1])
m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]
"""
# negbino 0.878
# zeroinfl 0.9909

### EPG - Marker haplotypes

# get the information of the statistics how to interpret the variation
m = vcat(repeat([1], size(epg,1)), repeat([2], size(epg,1)))
n = vcat(epg[:,:ratio], epg[:,:Allelfreq_Isr])
n = Int64.(floor.(n .* 100))
mn = DataFrame(hcat(m,n))

@rput mn
R"""
library(MASS)
library(pscl)
mn[,1] = as.factor(mn[,1])
m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]
"""
# negbino 0.683
# zeroinfl 0.951

### EPS - Single SNP level

# get the information of the statistics how to interpret the variation
m = vcat(repeat([1], size(eps,1)), repeat([2], size(eps,1)))
n = vcat(eps[:,:ratio], eps[:,:Allelfreq_Isr])
n = Int64.(floor.(n .* 100))
mn = DataFrame(hcat(m,n))

@rput mn
R"""
library(MASS)
library(pscl)
mn[,1] = as.factor(mn[,1])
m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]
"""
# negbino 0.0236
# zeroinfl 0.1687



######## casp to MACE comp with stats test #############################################

### EP - gene haplotypes

# get the information of the statistics how to interpret the variation
m = vcat(repeat([1], size(m1,1)), repeat([2], size(m1,1)))
n = vcat(m1[:,:ratio], m1[:,:Allelfreq_Isr])
n = Int64.(floor.(n .* 100))
mn = DataFrame(hcat(m,n))

@rput mn
R"""
library(MASS)
library(pscl)
mn[,1] = as.factor(mn[,1])
m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]
"""
# negbino 0.9020
# zeroinfl  0.6729

### EPG - Marker haplotypes

# get the information of the statistics how to interpret the variation
m = vcat(repeat([1], size(m3,1)), repeat([2], size(m3,1)))
n = vcat(m3[:,:ratio], m3[:,:Allelfreq_Isr])
n = Int64.(floor.(n .* 100))
mn = DataFrame(hcat(m,n))

@rput mn
R"""
library(MASS)
library(pscl)
mn[,1] = as.factor(mn[,1])
m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]
"""
# negbino 0.5605
# zeroinfl 0.9881

### EPS - Single SNP level

# get the information of the statistics how to interpret the variation
m = vcat(repeat([1], size(m2,1)), repeat([2], size(m2,1)))
n = vcat(m2[:,:ratio], m2[:,:Allelfreq_Isr])
n = Int64.(floor.(n .* 100))
mn = DataFrame(hcat(m,n))

@rput mn
R"""
library(MASS)
library(pscl)
mn[,1] = as.factor(mn[,1])
m1 =  summary(zeroinfl(x2~x1, mn, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]
"""
# negbino 0.67072
# zeroinfl 0.6894


###################################################################
## merge casp to mace to wgs in one table in the 3 levels
m1; ep;
names!(ep,Symbol.(string.("WGS_",String.(names(ep)))))
names!(m1,Symbol.(string.("MACE_",String.(names(m1)))))
mep = join(m1, ep, on= :MACE_Gene_ID => :WGS_Gene_ID )
names!(mep, vcat(names(eps)[1:6],names(mep)[7:18]))
delete!(mep, 13:17)
cor(mep[!,:WGS_Allelfreq_Isr], mep[!,:MACE_Allelfreq_Isr]) # 0.7652865005480041
CSV.write("CASP_to_WGS_to_MACE_physical_map.txt", mep;delim="\t")

m3; epg;
names!(epg,Symbol.(string.("WGS_",String.(names(epg)))))
names!(m3,Symbol.(string.("MACE_",String.(names(m3)))))
meg = join(m3, epg, on= :MACE_Gene_ID => :WGS_Gene_ID )
names!(meg, vcat(names(eps)[1:6],names(meg)[7:20]))
delete!(meg, 13:14)
cor(meg[!,:WGS_Allelfreq_Isr], meg[!,:MACE_Allelfreq_Isr]) # 0.6782170515970949
CSV.write("CASP_to_WGS_to_MACE_genetic_map.txt", meg; delim="\t")


m2; eps;
nam = names(eps)[1:6]
names!(eps,Symbol.(string.("WGS_",String.(names(eps)))))
names!(m2,Symbol.(string.("MACE_",String.(names(m2)))))
mes = join(m2, eps, on= :MACE_Gene_ID => :WGS_Gene_ID )
names!(mes, vcat(nam,names(mes)[7:22]))
delete!(mes, 18:28)
cor(mes[!,:WGS_Allelfreq_Isr], mes[!,:MACE_Allelfreq_Isr]) # 0.6351696572680174
CSV.write("CASP_to_WGS_to_MACE_single_SNPs.txt", mes; delim="\t")

##
############################################ correlate the relicates of the pools for each treatment ####
# WGS gene hap
pf = join(freq2[:F22P1KEP_Haplotypes], freq2[:F22P1K1_Haplotypes], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.8626669206759752 - 61495 haplotype genes

# test if they are equal
qwe = DataFrame(hcat(vcat(freq2[:F22P1KEP_Haplotypes].Allelfreq_Isr, freq2[:F22P1K1_Haplotypes].Allelfreq_Isr), vcat( repeat(1:1,size(freq2[:F22P1KEP_Haplotypes],1)),  repeat(2:2,size(freq2[:F22P1K1_Haplotypes],1)))))
qwe = qwe[.!(isnan.(qwe.x1)),:]
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100))
qwe[!,:x2] .= Int64.(qwe.x2)
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))
"""


# WGS marker hap
pf = join(freq[:F22P1KEP_Haplotypes], freq[:F22P1K1_Haplotypes], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.9589965371036292 - 5941 marker

# test if they are equal
qwe = DataFrame(hcat(vcat(freq[:F22P1KEP_Haplotypes].Allelfreq_Isr, freq[:F22P1K1_Haplotypes].Allelfreq_Isr), vcat( repeat(1:1,size(freq[:F22P1KEP_Haplotypes],1)),  repeat(2:2,size(freq[:F22P1K1_Haplotypes],1)))))
qwe = qwe[.!(isnan.(qwe.x1)),:]
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100))
qwe[!,:x2] .= Int64.(qwe.x2)
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))
"""


# WGS SNP level
cor(.!(isnan.(freq2[:F22P1KEP][!,:Allelfreq_Isr])), .!(isnan.(freq2[:F22P1K1][!,:Allelfreq_Isr]))) # 0.04096303794974753

# test if they are equal
qwe = DataFrame(hcat(vcat(freq[:F22P1KEP].Allelfreq_Isr, freq[:F22P1K1].Allelfreq_Isr), vcat( repeat(1:1,size(freq[:F22P1KEP],1)),  repeat(2:2,size(freq[:F22P1K1],1)))));
qwe = qwe[.!(isnan.(qwe.x1)),:];
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100));
qwe[!,:x2] .= Int64.(qwe.x2);
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))
"""


# WGS contig level
pf = join(freq[:F22P1KEP_Contigs], freq[:F22P1K1_Contigs], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.9855899118701991 - 485 contigs

# test if they are equal
qwe = DataFrame(hcat(vcat(freq[:F22P1KEP_Contigs].Allelfreq_Isr, freq[:F22P1K1_Contigs].Allelfreq_Isr), vcat( repeat(1:1,size(freq[:F22P1KEP_Contigs],1)),  repeat(2:2,size(freq[:F22P1K1_Contigs],1)))))
qwe = qwe[.!(isnan.(qwe.x1)),:]
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100))
qwe[!,:x2] .= Int64.(qwe.x2)
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))
"""





# MACE gene hap level
pf = join(fepm[:f22epkd_Haplotypes], fepm[:f22iiikd_Haplotypes], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.5681866898448314
pf = join(fepm[:f22epkd_Haplotypes], fepm[:f22iikd_Haplotypes], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.5892058387556539
pf = join(fepm[:f22iiikd_Haplotypes], fepm[:f22iikd_Haplotypes], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.6373534985200415

# test if they are equal
qwe = DataFrame(hcat(vcat(fepm[:f22iikd_Haplotypes].Allelfreq_Isr, fepm[:f22epkd_Haplotypes].Allelfreq_Isr, fepm[:f22iiikd_Haplotypes].Allelfreq_Isr), vcat( repeat(1:1,size(fepm[:f22epkd_Haplotypes],1)),  repeat(2:2,size(fepm[:f22iikd_Haplotypes],1)),  repeat(3:3,size(fepm[:f22iiikd_Haplotypes],1)))))
qwe = qwe[.!(isnan.(qwe.x1)),:]
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100))
qwe[!,:x2] .= Int64.(qwe.x2)
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))

"""


# MACE Marker hap level
pf = join(fegm[:f22epkd_Haplotypes], fegm[:f22iiikd_Haplotypes], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.6568517127413341
pf = join(fegm[:f22epkd_Haplotypes], fegm[:f22iikd_Haplotypes], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.6720309090400798
pf = join(fegm[:f22iiikd_Haplotypes], fegm[:f22iikd_Haplotypes], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.7060260585350422

# test if they are equal
qwe = DataFrame(hcat(vcat(fegm[:f22iikd_Haplotypes].Allelfreq_Isr, fegm[:f22epkd_Haplotypes].Allelfreq_Isr, fegm[:f22iiikd_Haplotypes].Allelfreq_Isr), vcat( repeat(1:1,size(fegm[:f22epkd_Haplotypes],1)),  repeat(2:2,size(fegm[:f22iikd_Haplotypes],1)),  repeat(3:3,size(fegm[:f22iiikd_Haplotypes],1)))))
qwe = qwe[.!(isnan.(qwe.x1)),:]
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100))
qwe[!,:x2] .= Int64.(qwe.x2)
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]

"""


# MACE Contig
pf = join(fepm[:f22epkd_Contigs], fepm[:f22iiikd_Contigs], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.9234367469696552
pf = join(fepm[:f22epkd_Contigs], fepm[:f22iikd_Contigs], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.9105693481412731
pf = join(fepm[:f22iiikd_Contigs], fepm[:f22iikd_Contigs], on=:Gene_ID, makeunique=true)
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.8432446256584535

# test if they are equal
qwe = DataFrame(hcat(vcat(fepm[:f22iikd_Contigs].Allelfreq_Isr, fepm[:f22epkd_Contigs].Allelfreq_Isr, fepm[:f22iiikd_Contigs].Allelfreq_Isr), vcat( repeat(1:1,size(fepm[:f22epkd_Contigs],1)),  repeat(2:2,size(fepm[:f22iikd_Contigs],1)),  repeat(3:3,size(fepm[:f22iiikd_Contigs],1)))))
qwe = qwe[.!(isnan.(qwe.x1)),:]
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100))
qwe[!,:x2] .= Int64.(qwe.x2)
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]

"""


# MACE SNP
cor(.!(isnan.(fepm[:f22epkd][!,:Allelfreq_Isr])), .!(isnan.(fepm[:f22iikd][!,:Allelfreq_Isr]))) # 0.49950203052128983
cor(.!(isnan.(fepm[:f22epkd][!,:Allelfreq_Isr])), .!(isnan.(fepm[:f22iiikd][!,:Allelfreq_Isr]))) # 0.49964210198020975
cor(.!(isnan.(fepm[:f22iiikd][!,:Allelfreq_Isr])), .!(isnan.(fepm[:f22iikd][!,:Allelfreq_Isr]))) # 0.5356893469124264

# test if they are equal
qwe = DataFrame(hcat(vcat(fepm[:f22iikd].Allelfreq_Isr, fepm[:f22epkd].Allelfreq_Isr, fepm[:f22iiikd].Allelfreq_Isr), vcat( repeat(1:1,size(fepm[:f22epkd],1)),  repeat(2:2,size(fepm[:f22iikd],1)),  repeat(3:3,size(fepm[:f22iiikd],1)))))
qwe = qwe[.!(isnan.(qwe.x1)),:]
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100))
qwe[!,:x2] .= Int64.(qwe.x2)
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]

"""


#########################
# GBS single SNP
fep[:f22epkd][!,:ID] .= string.(fep[:f22epkd].Chr,fep[:f22iikd].Pos)
fep[:f22iikd][!,:ID] .= string.(fep[:f22iikd].Chr,fep[:f22iikd].Pos)
fep[:f22iiikd][!,:ID] .= string.(fep[:f22iiikd].Chr,fep[:f22iikd].Pos)
pf = join(fep[:f22epkd], fep[:f22iiikd], on=:ID, makeunique=true);
pf = pf[.&(.!(isnan.(pf.Allelfreq_Isr)), .!(isnan.(pf.Allelfreq_Isr_1))),:];
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.55203049471351
pf = join(fep[:f22epkd], fep[:f22iikd], on=:ID, makeunique=true);
pf = pf[.&(.!(isnan.(pf.Allelfreq_Isr)), .!(isnan.(pf.Allelfreq_Isr_1))),:];
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.5439171451459219
pf = join(fep[:f22iiikd], fep[:f22iikd], on=:ID, makeunique=true);
pf = pf[.&(.!(isnan.(pf.Allelfreq_Isr)), .!(isnan.(pf.Allelfreq_Isr_1))),:];
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.5468133946541784

# test if they are equal
qwe = DataFrame(hcat(vcat(fep[:f22iikd].Allelfreq_Isr, fep[:f22epkd].Allelfreq_Isr, fep[:f22iiikd].Allelfreq_Isr), vcat( repeat(1:1,size(fep[:f22epkd],1)),  repeat(2:2,size(fep[:f22iikd],1)),  repeat(3:3,size(fep[:f22iiikd],1)))))
qwe = qwe[.!(isnan.(qwe.x1)),:]
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100))
qwe[!,:x2] .= Int64.(qwe.x2)
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]

"""

# GBS gene hap
pf = join(fep[:f22epkd_Haplotypes], fep[:f22iiikd_Haplotypes], on=:Gene_ID, makeunique=true);
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.6039900225769809
pf = join(fep[:f22epkd_Haplotypes], fep[:f22iikd_Haplotypes], on=:Gene_ID, makeunique=true);
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.6089527640472575
pf = join(fep[:f22iiikd_Haplotypes], fep[:f22iikd_Haplotypes], on=:Gene_ID, makeunique=true);
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.5978768797915586

# test if they are equal
qwe = DataFrame(hcat(vcat(fep[:f22iikd_Haplotypes].Allelfreq_Isr, fep[:f22epkd_Haplotypes].Allelfreq_Isr, fep[:f22iiikd_Haplotypes].Allelfreq_Isr), vcat( repeat(1:1,size(fep[:f22epkd_Haplotypes],1)),  repeat(2:2,size(fep[:f22iikd_Haplotypes],1)),  repeat(3:3,size(fep[:f22iiikd_Haplotypes],1)))))
qwe = qwe[.!(isnan.(qwe.x1)),:]
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100))
qwe[!,:x2] .= Int64.(qwe.x2)
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))
"""

# GBS contig correlation
pf = join(feg[:f22epkd_Contigs], feg[:f22iiikd_Contigs], on=:Gene_ID, makeunique=true);
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.9504327989176334
pf = join(feg[:f22epkd_Contigs], feg[:f22iikd_Contigs], on=:Gene_ID, makeunique=true);
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.9608036020851068
pf = join(feg[:f22iiikd_Contigs], feg[:f22iikd_Contigs], on=:Gene_ID, makeunique=true);
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.9588921334965985

# contigs equal?
qwe = DataFrame(hcat(vcat(feg[:f22iikd_Contigs].Allelfreq_Isr, feg[:f22epkd_Contigs].Allelfreq_Isr, feg[:f22iiikd_Contigs].Allelfreq_Isr), vcat( repeat(1:1,size(feg[:f22epkd_Contigs],1)),  repeat(2:2,size(feg[:f22iikd_Contigs],1)),  repeat(3:3,size(feg[:f22iiikd_Contigs],1)))))
qwe = qwe[.!(isnan.(qwe.x1)),:]
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100))
qwe[!,:x2] .= Int64.(qwe.x2)
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]

"""

# GBS Marker hap level
pf = join(feg[:f22epkd_Haplotypes], feg[:f22iiikd_Haplotypes], on=:Gene_ID, makeunique=true);
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.7568354441216598
pf = join(feg[:f22epkd_Haplotypes], feg[:f22iikd_Haplotypes], on=:Gene_ID, makeunique=true);
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.7443347217045672
pf = join(feg[:f22iiikd_Haplotypes], feg[:f22iikd_Haplotypes], on=:Gene_ID, makeunique=true);
cor(pf[!,:Allelfreq_Isr], pf[!,:Allelfreq_Isr_1]) # 0.7608423979985225

# marker haplotypes equal?
qwe = DataFrame(hcat(vcat(feg[:f22iikd_Haplotypes].Allelfreq_Isr, feg[:f22epkd_Haplotypes].Allelfreq_Isr, feg[:f22iiikd_Haplotypes].Allelfreq_Isr), vcat( repeat(1:1,size(feg[:f22epkd_Haplotypes],1)),  repeat(2:2,size(feg[:f22iikd_Haplotypes],1)),  repeat(3:3,size(feg[:f22iiikd_Haplotypes],1)))))
qwe = qwe[.!(isnan.(qwe.x1)),:]
qwe[!,:x1] .= Int64.(floor.(qwe.x1 * 100))
qwe[!,:x2] .= Int64.(qwe.x2)
@rput qwe
R"""
qwe$x2 = as.factor(qwe$x2)
library(MASS)
library(pscl)
m1 =  summary(zeroinfl(x1~x2, qwe, dist = "negbin"))
negbino = m1$coefficients$count[,4]
zeroinf = m1$coefficients$zero[,4]

"""


################################################################################################################
## test if the extension size has an impact on the accuracy of the freq estimation

# physical map
casp = CSV.read("CASP_to_WGS_to_MACE_physical_map.txt", delim="\t", header=true)
loc = CSV.read("Location.txt", delim="\t",header=true)
loc = unique(loc[!,3:12])
set = loc[indexin(casp[!,:Gene_ID], loc[!,:Gene_ID]),:]
casp[!,:extensionSize] .= set[!,:end] .- set[!,:start]
casp[!,:AbweichungWGS] .= casp[!,:ratio] .- casp[!,:WGS_Allelfreq_Isr]
casp[!,:AbweichungMACE] .= casp[!,:ratio] .- casp[!,:MACE_Allelfreq_Isr]
cor(casp[!,:extensionSize],casp[!,:AbweichungWGS]) # -0.179849237356019
cor(casp[!,:extensionSize],casp[!,:AbweichungMACE]) # 0.5727611492444719
cor(casp[!,:extensionSize],casp[!,:WGS_Readcount]) # 0.8075959165963394
cor(casp[!,:extensionSize],casp[!,:MACE_Readcount]) # -0.29144861378808595

 # calculate the regression of the
lm(@formula(AbweichungMACE~extensionSize+Pos), casp)
#                 Estimate   Std. Error   t value  Pr(>|t|)     Lower 95%    Upper 95%
#───────────────────────────────────────────────────────────────────────────────────────
#(Intercept)    -0.0375923    0.0177646    -2.11614    0.0504  -0.0752516    6.6912e-5
#extensionSize   2.24777e-7   7.18806e-8    3.12708    0.0065   7.23965e-8   3.77157e-7
#Pos             5.03813e-11  3.85392e-11   1.30727    0.2096  -3.13182e-11  1.32081e-10

lm(@formula(AbweichungWGS~extensionSize+Pos), casp)
#                 Estimate   Std. Error    t value  Pr(>|t|)     Lower 95%    Upper 95%
#────────────────────────────────────────────────────────────────────────────────────────
#(Intercept)     0.0188801    0.0198432     0.951464    0.3555  -0.0231857    0.0609459
#extensionSize  -6.67897e-8   8.02915e-8   -0.831841    0.4177  -2.37e-7      1.03421e-7
#Pos            -2.51521e-11  4.30487e-11  -0.58427     0.5672  -1.16411e-10  6.61072e-11

# genetic map
casp = CSV.read("CASP_to_WGS_to_MACE_genetic_map.txt", delim="\t", header=true)
loc = CSV.read("Marker_reference_all.txt", delim="\t",header=true)
set = loc[indexin(casp[!,:MACE_Gene_ID_1], loc[!,:Marker]),:]
casp[!,:extensionSize] .= set[!,:end] .- set[!,:start]
casp[!,:AbweichungWGS] .= casp[!,:ratio] .- casp[!,:WGS_Allelfreq_Isr]
casp[!,:AbweichungMACE] .= casp[!,:ratio] .- casp[!,:MACE_Allelfreq_Isr]
cor(casp[!,:extensionSize],casp[!,:AbweichungWGS]) # -0.252841973512165
cor(casp[!,:extensionSize],casp[!,:AbweichungMACE]) # -0.10210975909001506
cor(casp[!,:extensionSize],casp[!,:WGS_Readcount]) # 0.9733088227767234
cor(casp[!,:extensionSize],casp[!,:MACE_Readcount]) # -0.38725339614229487

lm(@formula(AbweichungMACE~extensionSize+Pos), casp)
#                  Estimate   Std. Error    t value  Pr(>|t|)     Lower 95%    Upper 95%
#────────────────────────────────────────────────────────────────────────────────────────
#(Intercept)    -0.050807     0.0286647    -1.77246     0.1267  -0.120947     0.019333
#extensionSize  -1.0453e-8    3.8706e-8    -0.270061    0.7962  -1.05163e-7   8.42572e-8
#Pos             1.05231e-10  6.16981e-11   1.70558     0.1390  -4.57388e-11  2.56201e-10

lm(@formula(AbweichungWGS~extensionSize+Pos), casp)
#                  Estimate   Std. Error    t value  Pr(>|t|)     Lower 95%    Upper 95%
#────────────────────────────────────────────────────────────────────────────────────────
#(Intercept)     0.0461012    0.033087      1.39333     0.2129  -0.0348597    0.127062
#extensionSize  -3.12856e-8   4.46774e-8   -0.700256    0.5100  -1.40607e-7   7.8036e-8
#Pos            -6.36513e-11  7.12166e-11  -0.89377     0.4059  -2.37912e-10  1.10609e-10


##
# get the count of snp and haplotypes for the casp comaprisen
# WGS
size(freq2[:Location],1) # snp
length(unique(freq2[:Location].Gene_ID)) # genes
length(unique(freq[:Location].Marker)) # marker
length(unique(freq[:Location].Contigs)) # contigs
# MACE
size(fepm[:Location],1) # snp
length(unique(fepm[:Location].Gene_ID)) # genes
length(unique(fegm[:Location].Marker)) # marker
length(unique(fepm[:Location].Contigs)) # contigs
#GBS
size(fep[:Location],1) # snp
length(unique(fep[:Location].Gene_ID)) # genes
length(unique(feg[:Location].Marker)) # marker
length(unique(feg[:Location].Contigs)) # contigs

# median read count
#WGS
median(vcat(freq[:F22P1KEP].Readcount, freq[:F22P1K1].Readcount)) # snp median over all replicates
median(vcat(freq2[:F22P1KEP_Haplotypes].Readcount, freq2[:F22P1K1_Haplotypes].Readcount)) # gene
median(vcat(freq[:F22P1KEP_Haplotypes].Readcount, freq[:F22P1K1_Haplotypes].Readcount)) # marker
median(vcat(freq[:F22P1KEP_Contigs].Readcount, freq[:F22P1K1_Contigs].Readcount)) # contig

#MACE
median(vcat(fepm[:f22iikd].Readcount, fepm[:f22iiikd].Readcount, fepm[:f22epkd].Readcount)) # snp median over all replicates
median(vcat(fepm[:f22iikd_Haplotypes].Readcount, fepm[:f22iiikd_Haplotypes].Readcount, fepm[:f22epkd_Haplotypes].Readcount)) # gene
median(vcat(fegm[:f22iikd_Haplotypes].Readcount, fegm[:f22iiikd_Haplotypes].Readcount, fegm[:f22epkd_Haplotypes].Readcount))
median(vcat(fepm[:f22iikd_Contigs].Readcount, fepm[:f22iiikd_Contigs].Readcount, fepm[:f22epkd_Contigs].Readcount))
#GBS
median(vcat(fep[:f22iikd].Readcount, fep[:f22iiikd].Readcount, fep[:f22epkd].Readcount)) # snp median over all replicates
median(vcat(fep[:f22iikd_Haplotypes].Readcount, fep[:f22iiikd_Haplotypes].Readcount, fep[:f22epkd_Haplotypes].Readcount)) # gene
median(vcat(feg[:f22iikd_Haplotypes].Readcount, feg[:f22iiikd_Haplotypes].Readcount, feg[:f22epkd_Haplotypes].Readcount))
median(vcat(feg[:f22iikd_Contigs].Readcount, feg[:f22iiikd_Contigs].Readcount, feg[:f22epkd_Contigs].Readcount))


## stat values over the whole genome

###################################################
##WGS
median(vcat(freq[:F22P1K1].Readcount, freq[:F22P1KEP].Readcount)) # snp
median(vcat(freq2[:F22P1K1_Haplotypes].Readcount, freq2[:F22P1KEP_Haplotypes].Readcount)) # gene
median(vcat(freq[:F22P1K1_Haplotypes].Readcount, freq[:F22P1KEP_Haplotypes].Readcount)) # marker
median(vcat(freq[:F22P1K1_Contigs].Readcount, freq[:F22P1KEP_Contigs].Readcount)) # contig

mean(vcat(freq[:F22P1K1].Readcount, freq[:F22P1KEP].Readcount)) # snp
mean(vcat(freq2[:F22P1K1_Haplotypes].Readcount, freq2[:F22P1KEP_Haplotypes].Readcount)) # gene
mean(vcat(freq[:F22P1K1_Haplotypes].Readcount, freq[:F22P1KEP_Haplotypes].Readcount)) # marker
mean(vcat(freq[:F22P1K1_Contigs].Readcount, freq[:F22P1KEP_Contigs].Readcount)) # contig

var(vcat(freq[:F22P1K1].Readcount, freq[:F22P1KEP].Readcount)) # snp
var(vcat(freq2[:F22P1K1_Haplotypes].Readcount, freq2[:F22P1KEP_Haplotypes].Readcount)) # gene
var(vcat(freq[:F22P1K1_Haplotypes].Readcount, freq[:F22P1KEP_Haplotypes].Readcount)) # marker
var(vcat(freq[:F22P1K1_Contigs].Readcount, freq[:F22P1KEP_Contigs].Readcount)) # contig

sw = vcat(freq[:F22P1KEP][.!(isnan.(freq[:F22P1KEP].Allelfreq_Isr)),[9,11]], freq[:F22P1K1][.!(isnan.(freq[:F22P1K1].Allelfreq_Isr)),[9,11]]); # snp
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.063827046190847

sw = vcat(freq2[:F22P1KEP_Haplotypes], freq2[:F22P1K1_Haplotypes]); # gene
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.06407028048709143
median(sw.SNPcount)
mean(sw.SNPcount)

sw = vcat(freq[:F22P1KEP_Haplotypes], freq[:F22P1K1_Haplotypes]); # marker
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.06397869848136238
median(sw.SNPcount)
mean(sw.SNPcount)

sw = vcat(freq[:F22P1KEP_Contigs], freq[:F22P1K1_Contigs]); # contig
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.063827046190847
median(sw.SNPcount)
mean(sw.SNPcount)

#############################
# GBS (feg / fep)
median(vcat(feg[:f22epkd].Readcount, feg[:f22iikd].Readcount, feg[:f22iiikd].Readcount)) # snp
median(vcat(fep[:f22epkd_Haplotypes].Readcount, fep[:f22iikd_Haplotypes].Readcount, fep[:f22iiikd_Haplotypes].Readcount )) # gene
median(vcat(feg[:f22epkd_Haplotypes].Readcount, feg[:f22iikd_Haplotypes].Readcount, feg[:f22iiikd_Haplotypes].Readcount)) # marker
median(vcat(feg[:f22epkd_Contigs].Readcount, feg[:f22iikd_Contigs].Readcount, feg[:f22iiikd_Contigs].Readcount)) # contig

mean(vcat(feg[:f22epkd].Readcount, feg[:f22iikd].Readcount, feg[:f22iiikd].Readcount)) # snp
mean(vcat(fep[:f22epkd_Haplotypes].Readcount, fep[:f22iikd_Haplotypes].Readcount, fep[:f22iiikd_Haplotypes].Readcount )) # gene
mean(vcat(feg[:f22epkd_Haplotypes].Readcount, feg[:f22iikd_Haplotypes].Readcount, feg[:f22iiikd_Haplotypes].Readcount)) # marker
mean(vcat(feg[:f22epkd_Contigs].Readcount, feg[:f22iikd_Contigs].Readcount, feg[:f22iiikd_Contigs].Readcount)) # contig

var(vcat(feg[:f22epkd].Readcount, feg[:f22iikd].Readcount, feg[:f22iiikd].Readcount)) # snp
var(vcat(fep[:f22epkd_Haplotypes].Readcount, fep[:f22iikd_Haplotypes].Readcount, fep[:f22iiikd_Haplotypes].Readcount )) # gene
var(vcat(feg[:f22epkd_Haplotypes].Readcount, feg[:f22iikd_Haplotypes].Readcount, feg[:f22iiikd_Haplotypes].Readcount)) # marker
var(vcat(feg[:f22epkd_Contigs].Readcount, feg[:f22iikd_Contigs].Readcount, feg[:f22iiikd_Contigs].Readcount)) # contig

sw = vcat(feg[:f22epkd][.!(isnan.(feg[:f22epkd].Allelfreq_Isr)),[9,11]], feg[:f22iikd][.!(isnan.(feg[:f22iikd].Allelfreq_Isr)),[9,11]], feg[:f22iiikd][.!(isnan.(feg[:f22iiikd].Allelfreq_Isr)),[9,11]]); # snp
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.063827046190847

sw = vcat(fep[:f22epkd_Haplotypes], fep[:f22iikd_Haplotypes], fep[:f22iiikd_Haplotypes]); # snp
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.06407028048709143
median(sw.SNPcount)
mean(sw.SNPcount)

sw = vcat(feg[:f22epkd_Haplotypes], feg[:f22iikd_Haplotypes], feg[:f22iiikd_Haplotypes]); # snp
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.06397869848136238
median(sw.SNPcount)
mean(sw.SNPcount)

sw = vcat(feg[:f22epkd_Contigs], feg[:f22iikd_Contigs], feg[:f22iiikd_Contigs]); # snp
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.063827046190847
median(sw.SNPcount)
mean(sw.SNPcount)


###########################
##MACE
median(vcat(fegm[:f22epkd].Readcount, fegm[:f22iikd].Readcount, fegm[:f22iiikd].Readcount)) # snp
median(vcat(fepm[:f22epkd_Haplotypes].Readcount, fepm[:f22iikd_Haplotypes].Readcount, fepm[:f22iiikd_Haplotypes].Readcount )) # gene
median(vcat(fegm[:f22epkd_Haplotypes].Readcount, fegm[:f22iikd_Haplotypes].Readcount, fegm[:f22iiikd_Haplotypes].Readcount)) # marker
median(vcat(fepm[:f22epkd_Contigs].Readcount, fepm[:f22iikd_Contigs].Readcount, fepm[:f22iiikd_Contigs].Readcount)) # contig

mean(vcat(fegm[:f22epkd].Readcount, feg[:f22iikd].Readcount, fegm[:f22iiikd].Readcount)) # snp
mean(vcat(fepm[:f22epkd_Haplotypes].Readcount, fepm[:f22iikd_Haplotypes].Readcount, fepm[:f22iiikd_Haplotypes].Readcount )) # gene
mean(vcat(fegm[:f22epkd_Haplotypes].Readcount, fegm[:f22iikd_Haplotypes].Readcount, fegm[:f22iiikd_Haplotypes].Readcount)) # marker
mean(vcat(fepm[:f22epkd_Contigs].Readcount, fepm[:f22iikd_Contigs].Readcount, fepm[:f22iiikd_Contigs].Readcount)) # contig

var(vcat(fegm[:f22epkd].Readcount, fegm[:f22iikd].Readcount, fegm[:f22iiikd].Readcount)) # snp
var(vcat(fepm[:f22epkd_Haplotypes].Readcount, fepm[:f22iikd_Haplotypes].Readcount, fepm[:f22iiikd_Haplotypes].Readcount )) # gene
var(vcat(fegm[:f22epkd_Haplotypes].Readcount, fegm[:f22iikd_Haplotypes].Readcount, fegm[:f22iiikd_Haplotypes].Readcount)) # marker
var(vcat(fepm[:f22epkd_Contigs].Readcount, fepm[:f22iikd_Contigs].Readcount, fepm[:f22iiikd_Contigs].Readcount)) # contig

sw = vcat(fegm[:f22epkd][.!(isnan.(fegm[:f22epkd].Allelfreq_Isr)),[9,11]], fegm[:f22iikd][.!(isnan.(fegm[:f22iikd].Allelfreq_Isr)),[9,11]], fegm[:f22iiikd][.!(isnan.(fegm[:f22iiikd].Allelfreq_Isr)),[9,11]]); # snp
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.063827046190847

sw = vcat(fepm[:f22epkd_Haplotypes], fepm[:f22iikd_Haplotypes], fepm[:f22iiikd_Haplotypes]); # snp
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.06407028048709143
median(sw.SNPcount)
mean(sw.SNPcount)

sw = vcat(fegm[:f22epkd_Haplotypes], fegm[:f22iikd_Haplotypes], fegm[:f22iiikd_Haplotypes]); # snp
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.06397869848136238
median(sw.SNPcount)
mean(sw.SNPcount)

sw = vcat(fepm[:f22epkd_Contigs], fepm[:f22iikd_Contigs], fepm[:f22iiikd_Contigs]); # snp
sum(sw.Readcount .* sw.Allelfreq_Isr) / sum(sw.Readcount) # 0.063827046190847
median(sw.SNPcount)
mean(sw.SNPcount)
