
combine_ids = function(v,sep=c('_', '.', '-', '&', '|',' ', ','))
{
  if(!inherits(v,'data.frame') && !inherits(v,'list'))
    return(as.character(v))
  
  v = lapply(v, as.character)
  ss=NULL
  while(is.null(ss))
  {
    for(s in sep)
    {
      if(!any(sapply(v, function(vec) any(grepl(s,vec,fixed=TRUE)))))
      {
        ss=s
        break
      }
    }
    sep = paste(rep(sep[1],2),collapse='')
  }
  v$sep=ss
  do.call(paste,v)
}


# prepare data from possibly different sources
# to do: also accept mst db
get_mml_data = function(dataSrc, qtpredicate, env, group)
{
  if(inherits(dataSrc, 'matrix'))
  {
    dat = dataSrc
    mode(dat) = 'integer'
    if(is.null(colnames(dat)))
      colnames(dat) = sprintf("item%06i",1:ncol(dat))

    if(!is.null(group) && length(group) != nrow(dat))
      stop(sprintf("Length of group (%i) is not equal to number of rows in data (%i)",
                   length(group),nrow(dat)))

    person_id = if(is.null(rownames(dat))) sprintf("p%09i",1:nrow(dat)) else rownames(dat)
    persons=tibble(person_id=person_id)
    if(!is.null(group))
      persons$group=group

  } else
  {
    dat = get_resp_matrix(dataSrc,qtpredicate,env)
    # to do: possibility to handle groups in resp_matrix should become part of dexter
    if(!is.null(group))
    {
      if(!is.character(group))
        stop("Group should be a character variable")

      if(inherits(dataSrc, "DBIConnection"))
      {
        # to do: allow booklet_id??
        g = dbGetQuery(dataSrc,sprintf("SELECT person_id, %s FROM dxPersons;",
                                       paste0(group, collapse=',')))
      } else if(inherits(dataSrc, 'data.frame'))
      {
        g = distinct(dataSrc[,c('person_id',group)], .data$person_id, .keep_all=TRUE)
      }
      g = g %>%
        mutate(person_id = factor(.data$person_id, levels=rownames(dat))) %>%
        filter(!is.na(.data$person_id)) %>%
        arrange(as.integer(.data$person_id))

      persons=g
      group = combine_ids(g[,group])
    } else
    {
      persons = tibble(person_id=rownames(dat))
    }
  }
  if(is.factor(group))
    group = droplevels(group)
  
  list(persons=persons,group=as.factor(group),dat=dat)
}


check_connected = function(design, fixed_items)
{
  flag = check_connected_c(design$items, design$groups, fixed_items)
  if(flag!=1)
  {
    if(flag==0)
      stop("Design is not connected and model is not identified.", call.=FALSE)
    # if msg is not overwritten the return code is unforeseen, should not happen
    msg = "design is not connected"
    if(flag %in% c(2,6))
    {
      msg = paste("Design is only connected by fixed items.",
                  "The calibration will be extremely sensitive to misfit of fixed items.")
    } else if(flag==4)
    {
      msg = paste("Design is only connected by assuming common population distributions.",
                  "The calibration will be extremely sensitive to sampling error.")

    }  else if(flag==8)
    {
      msg = paste("Design is only connected by fixed items in combination with assuming common population distributions.",
                  "The calibration will be extremely sensitive to misfit of fixed items AND sampling error.")
    }
    warning(paste("Unconnected design.",msg,"Ideally you should use a connected design."), call.=FALSE)
  }
}

