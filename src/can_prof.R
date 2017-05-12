

can_prof<-function( python_script,
         coordx,
         coordy,
         las_file,
         radius = 10,
         h.max = 80,
         h.min = 2,
         layer.thickness = 1,
         leaf_angle_dist = "spherical",
         max_return = 3, 
         method = "macarthur_horn", 
         out_filename = "temp.csv"
)
{
  
  system2('python', args=(sprintf('"%1$s" "%2$s" "%3$s" "%4$s" "%5$s" "%6$s" "%7$s" "%8$s" "%9$s" "%10$s" "%11$s" "%12$s"',
                                  python_script,
                                  coordx,
                                  coordy,
                                  las_file,
                                  radius,
                                  h.max,
                                  h.min,
                                  layer.thickness,
                                  leaf_angle_dist,
                                  max_return, 
                                  method, 
                                  out_filename)))
  return(NULL)
}

can_prof(
  python_script = "simple_canopy_profile_driver_command_line_function.py",
  coordx = 577601.951,
  coordy = 526741.5,
  las_file = "Carbon_plot_point_cloud_buffer.las",
  out_filename = "test.csv"
)