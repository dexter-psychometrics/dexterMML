
progress_width = function()
{
  if(!getOption("knitr.in.progress",FALSE) && getOption("dexter.progress", TRUE) && interactive())
  {
    getOption("width")
  } else {-1L}
}