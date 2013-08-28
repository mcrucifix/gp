define_mask <-
function() 
{
  blank <- function() matrix(NA, 64, 32 ) 

  tropforest <- blank ()
  tropforest [53:58,13:16] = 1
  tropforest [4:7,13:16] = 1

  ssahel <- blank()
  ssahel [(62:70) %% 64 + 1 ,19:21] = 1

  arctic <- blank()
  arctic = matrix(0,64,32)
  arctic [-(50:60), 31:32 ] = 1
  arctic [, 32 ] = 1
  arctic [28:42, 30:31 ] = 1
  arctic [1:10, 30:31 ] = 1
  arctic [62:64, 30:31 ] = 1

  necan <- blank()
  necan [40:47, 27:29] = 1

  siberia <- blank()
  siberia [12:22, 27:29] = 1
  siberia [12:25, 27:29] = 1
  siberia [26:31, 28:29] = 1
 
  rm('blank')

  return(as.list(environment()))
  }
