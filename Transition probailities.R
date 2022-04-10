# Transition probability calculator:
#functions
matrix_exp <- function(mat,n){
  if(!is.matrix(mat) | (nrow(mat)!=ncol(mat)) | n <0){
    return(0)
  }
  if(n == 0){
    return(diag(nrow(mat)))
  }
  matnew <- diag(nrow(mat))
  for(i in 1:n){
    matnew <- matnew %*% mat
  }
  matnew
}
# Transition probabilities
Transition_p = matrix (0, 2000, 4)
q =1
colnames (Transition_p) = c("wa" , "gl", "f0", "fA" )
# parameter space
for (gl in seq(3,15,1)){
  for(wa in seq (0.1,1,0.05)){
   # for (gr in seq(0,0.3,0.1)){
      
      #parameters
      # gl = 11
      # wa = 0.9
      gl_0 = gl-1
      gr = 0.0
      s_rounds = 20
      mu = 0
      mi = 0
      gen = 100
      
      
      gr_0 = 1
      mass = matrix (0, gl, gen+1)
      freq = matrix (0, gl, gen+1)
      allele_f = matrix (0, gl, 1)
      growth_rate = matrix (0, gl, 1)
      mean_allele_f = matrix (0, gen, 1)
      mean_growth = matrix (0, gen, 1)
      add_fit = matrix (0, gl, 1)
      drift_effect = matrix (0,gl,1)
      Growth_mat = matrix (0, gl, gl)
      selection_mat = matrix (0, gl, gl)
      drift_mut_mat = matrix (0, gl, gl)
      mig_mat = matrix (0, gl, gl)
      mass [gl_0,1] = 1
      freq [gl_0,1] = 100
      
      # calculating the frequencies and average fitness in each grain level
      
      allele_f = seq (0,1,1/(gl-1))
      allele_f_after_s = allele_f[]*wa/(allele_f[]*wa+(1-allele_f[]))
      mean_allele_f[1]=allele_f[gl_0] 
      # effect of group growth changes
      for (i in 1:gl) {
        Growth_mat [i,i] = gr_0 + (i-1)/(gl-1)*gr
        growth_rate[i] = gr_0 + (i-1)/(gl-1)*gr
      }
      #Growth_mat [1,1] = 1
      #Growth_mat [gl,gl] = 1
      # effect of individual frequency changes in each grain due to natural selection
      #version 1
      #for (i in 1:gl){
      #  add_fit [i] = (allele_f[i]*wa/(allele_f[i]*wa+1-allele_f[i]))-allele_f[i]
      #  j = add_fit[i]%/%(1/(gl-1))
      #  k = add_fit[i]%%(1/(gl-1))
      #  if (add_fit[i]<0){
      #    if (i>-j){
      #      selection_mat [i+j,i]= 1-(k/(1/(gl-1))) 
      #    }
      #    selection_mat [i+j+1,i]= k/(1/(gl-1))
      #  } else {
      #    if (i<(gl-j)){
      #      selection_mat [i+j+1,i]= k/(1/(gl-1))
      #    }
      #    selection_mat [i+j,i]= 1-(k/(1/(gl-1)))
      #  }
      #}
      #version_2 (birth_death model)
      for (j in 1:gl){
        selection_mat [j,j] = allele_f[j]*allele_f_after_s[j]+(1-allele_f[j])*(1- allele_f_after_s[j])
        if (j+1 <= gl){
          selection_mat [j+1,j] = (1-allele_f[j])*allele_f_after_s[j]
        }
        if (j-1 >= 1){
          selection_mat [j-1,j] = allele_f[j]*(1-allele_f_after_s[j])
        }
      }
      selection_mat_final = matrix_exp(selection_mat,s_rounds)
      # effect of drift_version_1
      #for (i in 2:(gl-1)){
      #  drift_effect [i] = allele_f[i]*(1-allele_f[i])/n
      #  drift_mut_mat [i,i]=1-2* drift_effect[i]
      #  drift_mut_mat [i+1,i] = drift_effect[i]
      #  drift_mut_mat [i-1,i] = drift_effect[i]
      #}
      # effect of drift _version_2
      for (i in 1:gl){
        for (j in 1:gl){
          drift_mut_mat [i,j] = allele_f[j]^(i-1)*(1-allele_f[j])^(gl-i)*choose(gl-1,i-1)
        }
      }
      # effect of drift _version_3 (binning. I decided to don't do it)
      #drift_mut_mat_n = matrix (0,n,n)
      #for (i in 1:n){
      #  for (j in 1:n){
      #    drift_mut_mat_n [i,j] = allele_f[j]^(i-1)*(1-allele_f[j])^(n-i+1)*choose(n,i-1)
      #  }
      #}
      #for (i in 1:gl){
      #  for (j in 1:gl){
      #    drift_mut_mat [i,j] = sum (drift_mut_mat_n[i:round(j+n/((gl-1))),j:round(i+n/((gl-1)))])
      #  }
      #}
      
      
      # effect of mutation (in both directions)
      drift_mut_mat = diag (gl)
      #the previous line is included because birth/death model has inherent drift
      drift_mut_mat [2,1] = mu
      drift_mut_mat [1,1] = 1-mu
      drift_mut_mat [gl-1,gl] = mu
      drift_mut_mat [gl,gl] = 1-mu
      
      
      # running the model for 'gen' generations + incorporating the effect of migration
      for (i in 2:gen){
        for (j in 1:gl){
          mig_mat [j,j] = allele_f[j]*mean_allele_f[i-1]+(1-allele_f[j])*(1-mean_allele_f[i-1])
          if (j+1 <= gl){
            mig_mat [j+1,j] = (1-allele_f[j])*mean_allele_f[i-1]
          }
          if (j-1 >= 1){
            mig_mat [j-1,j] = allele_f[j]*(1-mean_allele_f[i-1])
          }
        }
        mig_mat_final = diag (gl)
        for (k in 1:mi){
          mig_mat_final = mig_mat_final %*% mig_mat
        }
        mean_growth [i-1] = sum((freq [,i-1]/100)* growth_rate )
        freq [,i]= (Growth_mat/mean_growth[i-1]) %*% selection_mat_final %*% drift_mut_mat %*% (matrix_exp(mig_mat,mi)) %*% freq [,i-1]
        mean_allele_f[i] = sum (freq[,i]*allele_f[])/100
        
      }
      Transition_p [q,] = c(wa , gl,  freq [1,gen], freq[gl, gen] )
      q = q+1
    #}
  }
}

plot (Transition_p[,1],Transition_p[,4])
Transition_p <- as.data.frame(Transition_p)

Transition_p <- Transition_p %>% 
  mutate (Kimura_fA = (1 - exp(2 * (gl-1) * (wa - 1) * (gl-2)/(gl-1)))/(1 - exp(2 * (gl-1) * (wa - 1) ))) %>% 
  mutate (Kimura_fa = (1 - exp(2  * ((1 / wa) - 1))) / (1 - exp(2 * (gl-1) * ((1 / wa) - 1))))

library ("ggplot2")
ggplot (Transition_p, aes (x = wa , y = fA , color = gl))+
  geom_point()


