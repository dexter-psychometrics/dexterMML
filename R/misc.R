
progress_width = function()
{
  if(!getOption("knitr.in.progress",FALSE) && getOption("dexter.progress", FALSE) && interactive())
  {
    getOption("width")
  } else {-1L}
}