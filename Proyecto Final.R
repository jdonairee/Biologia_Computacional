# montserrat Navarrete 
library(seqinr)
library(ggplot2)
library(dplyr)

genes = read.fasta("sequenceWuhan.txt", forceDNAtolower = FALSE)
genesAust = read.fasta("5Australia(D.2).fasta", forceDNAtolower = FALSE)

datos = data.frame(
  muta = character(),
  cambioCodon = character(),
  cambioAmino = character(),
  pos = integer(),
  gen = character()
)
index = 1

tradCodon = c(
  "GAC" = "D", "GAU" = "D",
  "GAA" = "E", "GAG" = "E",
  "CGA" = "R", "CGC" = "R", "CGG" = "R", "CGU" = "R", "AGA" = "R", "AGG" = "R",
  "AAA" = "K", "AAG" = "K",
  "AAC" = "N", "AAU" = "N",
  "CAC" = "H", "CAU" = "H",
  "CAA" = "Q", "CAG" = "Q",
  "UCA" = "S", "UCC" = "S", "UCG" = "S", "UCU" = "S", "AGC" = "S", "AGU" = "S",
  "ACA" = "T", "ACC" = "T", "ACG" = "T", "ACU" = "T",
  "GCA" = "A", "GCC" = "A", "GCG" = "A", "GCU" = "A",
  "GGA" = "G", "GGC" = "G", "GGG" = "G", "GGU" = "G",
  "GUA" = "V", "GUC" = "V", "GUG" = "V", "GUU" = "V",
  "CCA" = "P", "CCC" = "P", "CCG" = "P", "CCU" = "P",
  "CUA" = "L", "CUC" = "L", "CUG" = "L", "CUU" = "L", "UUA" = "L", "UUG" = "L",
  "UUC" = "F", "UUU" = "F",
  "UAC" = "Y", "UAU" = "Y",
  "AUA" = "I", "AUC" = "I", "AUU" = "I",
  "AUG" = "M",
  "UGG" = "W",
  "UGC" = "C", "UGU" = "C",
  "-")

for (i in seq(1, length(genes))){
  
  if (i == 2) {
    next
  }
  
  geness = genes[[i]]
  attr1 = attr(geness, "Annot")
  vec = unlist(strsplit(attr1,"\\[|\\]|:|=|\\."));
  gene = vec[which(vec=="gene")+1]
  location = vec[which(vec=="location")+1]
  cat("\nGen:" , gene, location, "\n")
  
  secuencia = ""
  for(j in seq(1, length(genes[[i]]))){
    letra = genes[[i]][j]
    secuencia = paste(secuencia, letra, sep="")
  }
  
  rna = ""
  for(j in seq(1, nchar(secuencia))){
    letra = substr(secuencia, j, j)
    cambio = switch(letra, "T" = "U", letra)
    rna = paste(rna, cambio, sep="")
  }
  
  gen = genes[[i]]
  numNucleot = nchar(rna)
  numCodones = numNucleot/3
  cat("Wuhan   -->  Nucleótidos:", numNucleot, "   -->  Codones:", numCodones, "\n")
  
  for (k in seq(i, length(genesAust), 12)){
    
    secuenciaAust = ""
    for(j in seq(1, length(genesAust[[k]]))){
      letraAust = genesAust[[k]][j]
      secuenciaAust = paste(secuenciaAust, letraAust, sep="")
    }
    
    rnaAust = ""
    for(j in seq(1, nchar(secuenciaAust))){
      letraAust = substr(secuenciaAust, j, j)
      cambio = switch(letraAust, "T" = "U", letraAust)
      rnaAust = paste(rnaAust, cambio, sep="")
    }
    
    genAust = genesAust[[k]]
    numNucleotAust = nchar(rnaAust)
    numCodonesAust = numNucleotAust/3
    cat("Australia   -->  Nucleótidos:", numNucleotAust, "   -->  Codones:", numCodonesAust, "\n")
# ----------------------------------------------------------------
    
    rnaComp = rna
    rnaAustComp = rnaAust
    genComp = gen
    genAustComp = genAust
    
    if (nchar(rnaComp) != nchar(rnaAustComp)){ 
      
      rnaAli = unlist(strsplit(rnaComp, ""))
      rnaAliAust = unlist(strsplit(rnaAustComp, ""))
      
      rnaAli = c(" ", rnaAli)
      rnaAliAust = c(" ", rnaAliAust)
      
      m = matrix(data=0, nrow=length(rnaAliAust), ncol=length(rnaAli))
      m[1,] = seq(0,-2*(length(rnaAli)-1), -2)
      m[,1] = seq(0,-2*(length(rnaAliAust)-1), -2)
      for (fila in seq(2, length(rnaAliAust))){
        for (col in seq(2, length(rnaAli))){
          if (rnaAli[col]==rnaAliAust[fila]){
            d=m[fila-1, col-1] + 1
          } else {
            d=m[fila-1, col-1]-1
          }
          up = m[fila-1, col] -2
          left = m[fila, col-1] -2
          
          peso = max(d, up, left)
          m[fila, col] = peso
        }
      } 
      
      fila = length(rnaAliAust)
      col = length(rnaAli)
      newgen = c()
      newgenAust = c()
      while (fila > 1 || col > 1){
        if (fila > 1 && col > 1 && rnaAli[col] == rnaAliAust[fila]){
          newgen = c(rnaAli[col], newgen)
          newgenAust = c(rnaAliAust[fila], newgenAust)
          fila = fila-1
          col = col-1
        } 
        else{
          if (fila > 1 && col > 1){
            up = m[fila-1, col] - 2
            left = m[fila, col-1] - 2
            cambio = max(up, left)
          } else if (fila > 1){
            up = m[fila-1, col] - 2
            cambio = up
          } else if (col > 1){
            left = m[fila, col-1] - 2
            cambio = left
          }
          if (fila > 1 && cambio == up){
            newgen = c("-", newgen)
            newgenAust = c(rnaAliAust[fila], newgenAust)
            fila = fila-1
          }
          else if (col > 1 && cambio == left){
            newgen = c(rnaAli[col], newgen)
            newgenAust = c("-", newgenAust)
            col = col-1
          }
        }
      }
      rnaComp = paste(newgen, collapse="")
      rnaAustComp = paste(newgenAust, collapse="")
      
      genComp = unlist(strsplit(rnaComp, ""))
      genAustComp = unlist(strsplit(rnaAustComp, ""))
    }
    
# -----------------------------------------------------------------
    
    if (length(genComp) == length(genAustComp)){
      diff = which(genComp != genAustComp)
      if (length(diff) > 0){
        cat("  Se encontraron", length(diff), "diferencias \n")
        prevMutation = ""
        for (j in diff){
          muta = paste(genComp[j], " to ", genAustComp[j], sep="")
          inicio = j - ((j - 1) %% 3)
          codon1 = substr(rnaComp,inicio,inicio+2)
          codon2 = substr(rnaAustComp,inicio,inicio+2)
          cambio = paste(codon1, "to", codon2)
          numCodon = ((j - 1) %/% 3) + 1
          amino1 = tradCodon[codon1]
          amino2 = tradCodon[codon2]
          cambioAmino = paste(amino1, numCodon, amino2, sep="")
          if ((!is.na(amino1)) && (!is.na(amino2)) && (cambioAmino != prevMutation)){
            cat("  Mutación:   --> ", muta, " Codón:", codon1, "to", codon2, " Aminoacido:", cambioAmino, gene, numCodon, "\n")
            obs = list(muta, cambio, cambioAmino, j, gene)
            datos[index,] = obs
            index = index + 1
          }else {
            cat("  Diferencia con gap:   --> ", muta, " Codón:", codon1, "to", codon2,gene, numCodon, "\n")
          }
          prevMutation = cambioAmino
        }
      }
    }
  }
  
  aminoacidos = ""
  for(j in seq(1, numNucleot, 3)){
    codon = substr(rna, j, j+2)
    amino = tradCodon[codon]
    aminoacidos = paste(aminoacidos, amino, sep="")
  }
  aminoacidosAust = ""
  for(j in seq(1, numNucleot, 3)){
    codon = substr(rnaAust, j, j+2)
    amino = tradCodon[codon]
    aminoacidosAust = paste(aminoacidosAust, amino, sep="")
  }
  #cat("Wuhan\n   -->  ", aminoacidos, "\n")
  #cat("Australia\n   -->  ", aminoacidosAust, "\n")
}

# ==============================

dfMuta = filter(
  summarise(
    group_by(datos, muta),
    cuenta = n()),
  cuenta > 0.20 * sum(cuenta)
)
dfMuta = as.data.frame(dfMuta)
str(dfMuta)

p1 = ggplot(dfMuta)
p1 = p1 + aes(x=muta, y=cuenta, fill=muta, label=cuenta)
p1 = p1 + ggtitle("Cambio de Nucleótido")
p1 = p1 + labs(x="Mutation", y="Frecuencia", fill="cambio")
p1 = p1 + geom_bar(stat="identity")
p1 = p1 + geom_text(vjust=0)
p1

# ==============================

dfAmino = group_by(datos, gen, cambioAmino)

dfAmino = summarise(
  dfAmino,
  mutacion = first(muta),
  cambioCodon = first(cambioCodon),
  pos = first(pos),
  cuenta = n()
)

dfAmino = filter(
  dfAmino,
  cuenta > 0.20 * (length(genesAust)/12)
)

dfAmino = as.data.frame(dfAmino)
str(dfAmino)

p4 = ggplot(dfAmino)
p4 = p4 + aes(x=cambioAmino, y=cuenta, fill=cambioAmino, label=cuenta)
p4 = p4 + ggtitle("Cambio de Aminoácidos")
p4 = p4 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
p4 = p4 + geom_bar(stat="identity")
p4 = p4 + geom_text(vjust=0)
p4 = p4 + facet_grid(~gen, scales="free", space="free_x")
p4
