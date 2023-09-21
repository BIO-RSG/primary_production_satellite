library(DiagrammeR)
library(DiagrammeRsvg)
library(dplyr)
library(rsvg)

# DIAGRAM FOR PP AT TIME T
p1 <- DiagrammeR::grViz("
  digraph boxes_and_circles {

  # node definitions with substituted label text
  node [shape = rectangle]
  i [label = 'Solar zenith angle']
  h [label = 'Chl-a at depth z']
  e [label = 'PI parameters\nat depth z']
  j [label = 'Absorption and backscattering\nat depth z']
  k [label = 'I(z=0,λ)']
  {rank = same; h; e; k}
  l [label = 'Attenuation coefficient\nat depth z']
  {rank = same; l}
  m [label = 'I(z,λ)']
  {rank = same; m}
  n [label = 'P(z)']
  {rank = same; n}
  p [label = 'In situ PI parameters\nin water column']
  q [label = 'Satellite PAR']
  r [label = 'Chl-a profile in\nwater column']
  b [label = 'Latitude']
  c [label = 'Day of year']
  d [label = 'Time']
  {rank = same; p; q; r; b; c; d}
  a [label = 'PAR(z)']

  # connections
  b -> i
  c -> i
  d -> i
  p -> e [label = 'Linear interpolation\nto fill missing depths']
  e -> n
  q -> k [label = 'Scaling\nfactor']
  h -> j
  h -> n
  i -> k [label = 'Gregg-Carder\nlight model']
  j -> l
  k -> m
  l -> m
  a -> n
  m -> a [label = 'Integration over wavelength']
  r -> h [label = 'Linear interpolation\nto fill missing depths']

  }
  ", height = 1000)
p1 %>% export_svg %>% charToRaw %>% rsvg_png("images/pp_biov2_flowchart.png")



# # DIAGRAM FOR PP AT TIME T
# p1 <- DiagrammeR::grViz("
#   digraph boxes_and_circles {
# 
#   # node definitions with substituted label text
#   node [shape = rectangle]
#   i [label = 'Solar zenith angle']
#   h [label = 'Chl-a at depth z']
#   e [label = 'PI parameters\nat depth z']
#   j [label = 'Absorption and backscattering\nat depth z']
#   k [label = 'I(z=0,λ)']
#   {rank = same; h; e; k}
#   l [label = 'Attenuation coefficient\nat depth z']
#   {rank = same; l}
#   m [label = 'I(z,λ)']
#   {rank = same; m}
#   n [label = 'P(z)']
#   {rank = same; n}
#   o [label = 'Total primary production\nin the water column']
#   {rank = same; o}
#   p [label = 'In situ PI parameters\nat 2 depths']
#   q [label = 'Satellite PAR']
#   r [label = 'Chl-a profile in\nwater column']
#   b [label = 'Latitude']
#   c [label = 'Day of year']
#   d [label = 'Time']
#   {rank = same; p; q; r; b; c; d}
#   a [label = 'PAR(z)']
# 
#   # connections
#   b -> i
#   c -> i
#   d -> i
#   p -> e [label = 'Linear interpolation\nto fill missing depths']
#   e -> n
#   q -> k [label = 'Scaling\nfactor']
#   h -> j
#   h -> n
#   i -> k [label = 'Gregg-Carder\nlight model']
#   j -> l
#   k -> m
#   l -> m
#   a -> n
#   m -> a [label = 'Integration over wavelength']
#   n -> o [label = 'Integration over depth']
#   r -> h [label = 'Linear interpolation\nto fill missing depths']
# 
#   }
#   ", height = 1000)
# p1 %>% export_svg %>% charToRaw %>% rsvg_png("images/pp_biov2_flowchart.png")
