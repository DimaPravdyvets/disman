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
#function for N numbers
install.packages('rlist')
library(rlist)


#input should be in this format: (list(N11,N12,N13,N1N,N21,N22,N23,N2N,ETC),#of alleles)
#remembre that 12 and 21 is the same, so in case you need to, just divide those combinations by 2
# Constructing Quadratic Formula
result <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
    x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
    result = max(x_1,x_2)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    x = -b/(2*a)
  }
  else {"There are no real roots."} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
  b^2-4*a*c
}

a <- result(1,-2,1); 
a

N=((n(*n+1))/2)
n*n+n-2N=0


Ngen_greq <- function(L,n){
  N=0
  for (i in L){
    N=N+i
  }
  #print(floor(result(1,1,-2*length(L))))
  falafel = vector("list")
  step = 1
  while (step<=length(L)){
    p = (L[[step]]/N)
    falafel=c(falafel,sqrt(p))
    step=step+n+1
  }
  freq = vector("list")
  for (j in L){
    w = (j/N)
    freq=c(freq,w)
  }
  #print(do.call(sum,freq))
  #print(do.call(sum,fallel))
  print("frequencies genotypes:")
  print(paste(freq))
  print(paste("Sum of genotype frequencies: ",Reduce('+',freq)))
  
  print("frequencies alleles:")
  print(paste(falafel))
  print(paste("sum of alleles frequencies:",Reduce('+',falafel)))
}
Ngen_greq(list(1787,3037/2,3037/2,1305),2)
Ngen_greq(list(1787,3037/2,1305/2,1234/2,2343,3325/2,1567,2235/2,3321),3)
Ngen_greq(list(1787,3037/2,1305/2,1234/2,2343/2,3325,1567/2,2235/2,3321/2,3254,1458/2,2541/2,4444/2,4565/2,1211/2,1111),4)
Ngen_greq(list(11,12/2,13/2,14/2,15/2,21/2,22,23/2,24/2,25/2,31/2,32/2,33,34/2,35/2,41/2,42/2,43/2,44,45/2,51/2,52/2,53/2,54/2,55),5)




#MENDELIAN POPULATION
#RANDOM MATING -HWE
#One-gene two-alleles (X1 and X2). Genotypes: X11, X12 and X22
random_mating<-function(f11,f12,f22){
  p=f11+f12/2;
  q=1-p
  #mating pairs 
  f11_11=f11*f11
  f11_12=0#2*f11*f12
  f11_22=0#2*f11*f22
  f12_12=f12*f12
  f12_22=0#2*f12*f22
  f22_22=f22*f22
  #Genotype probabilities of progeny (offspring) f11_, f12_, f22_
  f11_=f11_11+f11_12/2+f12_12/4
  f12_=f11_12/2+f11_22+f12_12/2+f12_22/2
  f22_=f22_22+f12_22/2+f12_12/4
  #allele probabilitty
  p_=f11_+f12_/2
  q_=1-p_
  
  
  
  
  
  # output
  print(paste("Initial:", "P",p,"Q",q))
  print(paste("f11,f12,f22",f11_t,f12_t,f22_t))
  print(paste("Next generation:", "P",p_,"Q",q_))
  #print(f11,f12,f22)
  #print(p_,q_)
  #print(f11_,f12_,f22_)
  #print(p_*p_,2*p_*q_,q_*q_)
}
random_mating(0.25,0.5,0.25)

#Ex on this programm
#ex1

random_mating1<-function(f11,f12,f22){
  p=f11+f12/2;
  q=1-p
  #mating pairs 
  f11_11=f11*f11/sum(f11*f11,f12*f12,f22*f22)
  f12_12=f12*f12/sum(f11*f11,f12*f12,f22*f22)
  f22_22=f22*f12/sum(f11*f11,f12*f12,f22*f22)
  f11_12=0
  f11_22=0
  f12_22=0
  #Genotype probabilities of progeny (offspring) f11_, f12_, f22_
  f11_=f11_11+f11_12/2+f12_12/4
  f12_=f11_12/2+f11_22+f12_12/2+f12_22/2
  f22_=f22_22+f12_22/2+f12_12/4
  #allele probabilitty
  p_=f11_+f12_/2
  q_=1-p_
  
  
  
  
  
  # output
  #print(paste("Initial:", "P",p,"Q",q))
  #print(paste("f11,f12,f22",f11_,f12_,f22_))
  #print(paste("Next generation:", "P",p_,"Q",q_))
  #print(f11,f12,f22)
  #print(p_,q_)
  #print(f11_,f12_,f22_)
  #print(p_*p_,2*p_*q_,q_*q_)
  print(paste(f11_,f12_,f22_))
}
random_mating1(0.25,0.5,0.25)
random_mating1(0.333333333333333,0.333333333333333,0.5)



#segregation execrice (next slide)
