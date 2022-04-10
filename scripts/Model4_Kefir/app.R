#Javad's Shiny App
#packages
library(shiny)
library(tidyverse)
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
# Define UI for application that Model_4 which shows spread of a mutation in kefir community. 
ui <- fluidPage(
  headerPanel = headerPanel("BIOL336: Haploid model of selection"),
  # Application title
  titlePanel("Model_4"),
  # Sidebar with a slider input for number of bins 
  
  sidebarLayout(
    sidebarPanel(
      
      HTML("<p style='font-size:14px'><B>Frequency of allele A over time (p).</B>"),
      
      numericInput(inputId = "gl", label = "Grain levels", value = 101, 
                   min = 2, max = 200, step = 1),
      numericInput(inputId = "wa", label = "wa", value = 1, 
                   min = 0, max = 2, step = 0.01),
      numericInput(inputId = "gl_0", label = "Initial state (grain level)", value = 2, 
                   min = 0, max = 200, step = 1),
      numericInput(inputId = "s_rounds", label = "number of birth/death events per group selection", value = 2, 
                   min = 0, max = 200, step = 1),
      numericInput(inputId = "gr", label = "Max grain growth rate change", value = 0, 
                   min =-2 , max = 2, step = 0.01),
      numericInput(inputId = "mu", label = "mutation rate", value = 0.001, 
                   min = 0, max = 0.01, step = 0.0001),
      numericInput(inputId = "mi", label = "migration rate", value = 1, 
                   min = 0, max = 20, step = 1),
      numericInput(inputId = "gen", label = "number of generations", value = 100, 
                   min = 2, max = 1000, step = 10)
      
    ),
    
    mainPanel =  mainPanel(
      plotOutput(outputId = 'viz')
    )
  )
)
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$viz <- renderPlot({
    #parameters
    gl = input$gl
    wa = input$wa
    gl_0 = input$gl_0 
    s_rounds = input$s_rounds
    gr = input$gr
    mu = input$mu 
    mi = input$mi
    gen = input$gen
    #defining needed variables
    gr_0 = 1
    mass = matrix (0, gl, gen+1)
    freq = matrix (0, gl, gen+1)
    allele_f = matrix (0, gl, 1)
    allele_f_selection = matrix (0, gl, 1)
    mean_allele_f = matrix (0, gen, 1)
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
    allele_f_after_s = allele_f*wa/(allele_f*wa+(1-allele_f))
    mean_allele_f[1]=allele_f[gl_0] 
    # effect of group growth changes
    for (i in 1:gl) {
      Growth_mat [i,i] = gr_0 + (i-1)/(gl-1)*gr
    }
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
    #effect of drift _version_2
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
      mass [,i]= Growth_mat %*% selection_mat_final %*% drift_mut_mat %*% (matrix_exp(mig_mat,mi)) %*% mass [,i-1]
      freq [,i]= 100* mass [,i]/ sum(mass[,i])
      mean_allele_f[i] = sum (freq[,i]*allele_f[])/100
      
    }
    plot (allele_f,freq[,gen])
  })
} 
# Run the application 
shinyApp(ui = ui, server = server)