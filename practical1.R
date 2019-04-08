#function to estimate q and p of genotypes

gene_freq<-function(N11,N12,N22){
  N=N11+N12+N22;#N total sample
  p=(N11+N12/2)/N; #p and q allele frequency X1 and X2 respectively
  q=1-p;
  f11=N11/N;
  f12=N11/N;
  f22=N22/N;
  print(paste(p,"=p",q,'=q',f11,'=f11',f12,'=f12',f22,'=f22'));
  }
gene_freq(1787,3037,1305)

#MENDELIAN POPULATION
#RANDOM MATING -HWE
#One-gene two-alleles (X1 and X2). Genotypes: X11, X12 and X22
random_mating<-function(f11,f12,f22){
  p=f11+f12/2;
  q=1-p
  #mating pairs 
  f11_11=f11*f11
  f11_12=2*f11*f12
  f11_22=2*f11*f22
  f12_12=f12*f12
  f12_22=2*f12*f22
  f22_22=f22*f22
  #Genotype probabilities of progeny (offspring) f11_, f12_, f22_
  f11_=f11_11+f11_12/2+f12_12/4
  f12_=f11_12/2+f11_22+f12_12/2+f12_22/2
  f22_=f22_22+f12_22/2+f12_12/4
  #allele probabilitty
  p_=f11_+f12_/2
  q_=1-p_
  # output
  p
  q
  #print(f11,f12,f22)
  #print(p_,q_)
  #print(f11_,f12_,f22_)
  #print(p_*p_,2*p_*q_,q_*q_)
}
random_mating(0.5,0,0.5)
