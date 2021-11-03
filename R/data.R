

mml_pre = function(dataSrc, qtpredicate, env, group=NULL, sorted=TRUE)
{
  if(inherits(dataSrc, 'matrix'))
  {
    if(!is.null(qtpredicate)) stop("predicates are not supported for matrix dataSrc")
    
    dat = dataSrc
    mode(dat) = 'integer'
    if(is.null(colnames(dat)))
      colnames(dat) = sprintf("item%06i",1:ncol(dat))
    
    if(!is.null(group) && length(group) != nrow(dat))
      stop(sprintf("Length of group (%i) is not equal to number of rows in data (%i)",
                   length(group),nrow(dat)))
    
    person_id = if(is.null(rownames(dat))) 1:nrow(dat) else rownames(dat)
    if(is.null(group))
    {
      persons = tibble(person_id=person_id, c_group_nbr=0L)
      subgroups = tibble(c_group_nbr=0L, group_id='population', group_n = nrow(dat))
    } else
    {
      persons = tibble(person_id=person_id, group_id=factor(group), c_group_nbr = as.integer(.data$group_id)-1L)
      subgroups = count(persons, .data$c_group_nbr, name='group_n') %>%
        arrange(.data$c_group_nbr)
      subgroups$group_id = levels(persons$group_id)
    }

    rg = range(dat)
    if(rg[1] < 0L)
      stop('negative scores are not allowed')

    pre = mat_pre(dataSrc, rg[2])
    item_id = colnames(dat)
        
  } else
  {
    # aiaiaiai, an internal function
    # has to be done differently before we go to CRAN
    # to do: adapt dexter to export a get_responses_ with factors
    dat = dexter:::get_responses_(dataSrc,qtpredicate=qtpredicate,env=env, 
                                  columns=c('person_id','item_id','item_score',group))
    
    dat$person_id = factor(dat$person_id)
    dat$item_id = factor(dat$item_id)
    dat$item_score = as.integer(dat$item_score)
    
    if(!is.null(group))
    {
      if(!is.character(group))
        stop("Group should be a character vector")
      
      persons = distinct(dataSrc[,c('person_id',group)], .data$person_id, .keep_all=TRUE) 
      
      subgroups = persons[,group] %>% count(across(), name='group_n')
      subgroups$c_group_nbr = 0:(nrow(subgroups)-1L)
      
      persons = inner_join(persons, select(subgroups,-.data$group_n), by=group) %>%
        arrange(.data$person_id)

    } else
    {
      persons = tibble(person_id=levels(dat$person_id), c_group_nbr=0L)
      subgroups = tibble(c_group_nbr=0L, group_id='population', group_n = nlevels(dat$person_id))
    }

    rg = range(dat$item_score)
    if(anyNA(rg)) stop('item_score should not contain NA values')
    if(rg[1]<0L) stop('negative scores are not allowed')
    
    pre = df_pre(dat$person_id, dat$item_id, dat$item_score, rg[2], nlevels(dat$person_id), nlevels(dat$item_id), sorted)
    item_id = levels(dat$item_id)
    
  }
      
  
  list(item_id=item_id, persons=persons, groups=arrange(subgroups,.data$c_group_nbr), pre=pre)
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

