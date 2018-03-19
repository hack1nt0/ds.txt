setClass("dtm", slots = c(dl = "integer", df = "integer"), contains = "dgCMatrix")

setMethod("[", signature(x = "dtm", i = "index", j = "index", drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              super <- as(x, 'dgCMatrix')[i,j,drop = F]
              if (drop)
                  super
              else {
                  ans = as(super, 'dtm')
                  .Call('C_update_dtm_subset', PACKAGE = 'ds.txt', ans)
                  ans
              }
          })
setMethod("[", signature(x = "dtm", i = "index", j = "missing", drop = "logical"),
          function(x, i, j, ..., drop = TRUE) { #todo
              super <- as(x, 'dgCMatrix')[i,,drop = F]
              if (drop)
                  super
              else {
                  ans = as(super, 'dtm')
                  .Call('C_update_dtm_subset', PACKAGE = 'ds.txt', ans)
                  ans
              }
          })
setMethod("[", signature(x = "dtm", i = "missing", j = "index", drop = "logical"),
          function(x, i, j, ..., drop = TRUE) {
              super <- as(x, 'dgCMatrix')[,j,drop = F]
              if (drop)
                  super
              else {
                  ans = as(super, 'dtm')
                  .Call('C_update_dtm_subset', PACKAGE = 'ds.txt', ans)
                  ans
              }
          })
