DoEstRare=function(pheno, geno, position, genome.size, perm=NULL, alpha=NULL, c=NULL, autosomal=TRUE, gender=NULL){

  #check
  #phenotype
  if(!is.numeric(pheno)){
    stop("argument \"pheno\" requires a numeric vector")
  }
  tab=as.data.frame(table(pheno))
  if(nrow(tab)==2){
  }else{
    stop("The phenotype is not a binary trait")
  }

  if(length(pheno)!=nrow(geno)){
    stop("the \"pheno\" vector length is not equal to the number of rows in the \"geno\" matrix")
  }

  #position
  if(!is.numeric(position)){
    stop("argument \"position argument\" requires a numeric vector")
  }

  if(length(position)!=ncol(geno)){
    stop("the \"position\" vector length is not equal to the number of columns in the \"geno\" matrix")
  }

  #genome.size
  if(!is.numeric(genome.size)){
    stop("argument \"genome.size\" requires a positive numeric value")
  }else{
    if(genome.size<0){
      stop("argument \"genome.size\" requires a positive numeric value")
    }
  }
  #check if genome.size is correct

  #permutations
  if(is.null(perm)){
    if(is.null(alpha)| is.null(c)){
      stop("Please use arguments \"perm\" or \"alpha\" + \"c\"  for permutations")
    }
  }


  if(length(unique(position))<=1){
    stat_obs=NA
    p.value=NA
    warning("Only one rare variant position, NA value returned")
  }else{

    P=ncol(geno)   #number of variants
    N=nrow(geno)   #number of individuals
    from=1
    to=genome.size
    bw=bw.nrd0(position)

    stat_obs=.C("DoEstRare_stat",
                as.numeric(0),
                as.numeric(c(geno)),
                as.numeric(pheno),
                as.numeric(N),
                as.numeric(P),
                as.numeric(autosomal),
                as.numeric(gender),
                as.numeric(position),
                as.numeric(from),
                as.numeric(to) ,
                as.numeric(bw))[[1]]

    if(!is.null(perm)){
      stat_perm=rep(NA,perm)
      for(i in 1:perm){
        pheno_perm=sample(pheno)
        stat_perm[i]=.C("DoEstRare_stat",
                        as.numeric(0),
                        as.numeric(c(geno)),
                        as.numeric(pheno_perm),
                        as.numeric(N),
                        as.numeric(P),
                        as.numeric(autosomal),
                        as.numeric(gender),
                        as.numeric(position),
                        as.numeric(from),
                        as.numeric(to) ,
                        as.numeric(bw))[[1]]
      }
      p.value=(length(which(stat_perm>=stat_obs))+1)/(perm+1)

    }else{
      b=choose_b(alpha, c)
      r=choose_r(alpha, c)

      j=1
      Rij=0
      while(Rij<r & j<b){
        stat_perm=NA
        pheno_perm=sample(pheno)

        stat_perm=.C("DoEstRare_stat",
                     as.numeric(0),
                     as.numeric(c(geno)),
                     as.numeric(pheno_perm),
                     as.numeric(N),
                     as.numeric(P),
                     as.numeric(autosomal),
                     as.numeric(gender),
                     as.numeric(position),
                     as.numeric(from),
                     as.numeric(to) ,
                     as.numeric(bw))[[1]]

        if(stat_perm>=stat_obs){
          Rij=Rij+1
        }

        j=j+1
      }

      if(j<b){
        p.value=Rij/(j-1)
      }else {
        p.value=(Rij+1)/(b+1)
      }
    }
  }


  res_DoEstRare=list(stat=stat_obs, p.value=p.value)
  return(res_DoEstRare)
}
