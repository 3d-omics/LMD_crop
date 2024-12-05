# Parameters

cryosection = "G007bI105A"
overview_image="input/241113G007bI105post.jpg"

cryosection = "G006bI203A"
overview_image="input/241113G006bI203post.jpg"

cryosection = "G082bI105A"
overview_image="input/241104G082bI105post.jpg"

cryosection = "G082bI105B"
overview_image="input/241104G082bI105post.jpg"

# Crop width and height
width=1000
height=1000

# Membrane dimensions in microns
membrane_micron_width=16000
membrane_micron_height=45000

membrane_pixel_width=2180
membrane_pixel_height=5780

# Load libraries

suppressPackageStartupMessages(library(striprtf))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(rairtable))
suppressPackageStartupMessages(library(magick))

# Declare slide position coordinates

slide1_tl <- c(250,2212) # (x, y) for top-left
slide1_br <- c(25660,77470)    # (x, y) for bottom-right
membrane1_tl <- c(4962,8698) # (x, y) for top-left
membrane1_br <- c(21013,53732) # (x, y) for bottom-right

slide2_tl <- c(28983,2004) # (x, y) for top-left
slide2_br <- c(54565,77541) # (x, y) for bottom-right
membrane2_tl <- c(33689,8698) # (x, y) for top-left
membrane2_br <- c(49650,53658) # (x, y) for bottom-right

slide3_tl <- c(64978,1938) # (x, y) for top-left
slide3_br <- c(90644,77372) # (x, y) for bottom-right
membrane3_tl <- c(69768,8658) # (x, y) for top-left
membrane3_br <- c(85819,53692) # (x, y) for bottom-right

slide4_tl <- c(93780,1946) # (x, y) for top-left
slide4_br <- c(119690,77465) # (x, y) for bottom-right
membrane4_tl <- c(98718,8682) # (x, y) for top-left
membrane4_br <- c(114760,53772) # (x, y) for bottom-right

# Collect input coordinates from airtable (requires personal token)
input_data <- airtable("4-MSE-Info", "appKakM1bnKSekwuW") %>% #get base ID from Airtable browser URL
  read_airtable(., fields = c("ID","Xcoord","Ycoord","size","shape","cryosection_text","SampleType"), id_to_col = TRUE) %>% #get 3 columns from MAGs table
  filter(cryosection_text == cryosection) %>% 
  select(-airtable_record_id) %>% 
  as_tibble() %>% 
  unnest(cols = c(Xcoord, Ycoord, size, shape)) %>%
  rename(microsample=ID,cryosection=cryosection_text)




#Infer slide position from coordinates and declare references
microsample_centroid <- input_data %>% 
  summarise(mean_x=mean(Xcoord),
            mean_y=mean(Ycoord))

if(microsample_centroid$mean_x > membrane1_tl[1] && microsample_centroid$mean_x < membrane1_br[1]){
  slide_position=1
  membrane_tl=membrane1_tl
  membrane_br=membrane1_br
}else if(microsample_centroid$mean_x > membrane2_tl[1] && microsample_centroid$mean_x < membrane2_br[1]){
  slide_position=2
  membrane_tl=membrane2_tl
  membrane_br=membrane2_br
}else if(microsample_centroid$mean_x > membrane3_tl[1] && microsample_centroid$mean_x < membrane3_br[1]){
  slide_position=3
  membrane_tl=membrane3_tl
  membrane_br=membrane3_br
}else if(microsample_centroid$mean_x > membrane4_tl[1] && microsample_centroid$mean_x < membrane4_br[1]){
  slide_position=4
  membrane_tl=membrane4_tl
  membrane_br=membrane4_br
}else{
  slide_position=NA
}

# Calculate resolution
resolution_x=membrane_micron_width / membrane_pixel_width
resolution_x=7.4 # this line added to "correct" inferred resolution
resolution_y=membrane_micron_height / membrane_pixel_height

# Calculate microsample coordinates in overview image
input_data <- input_data %>% 
  mutate(Xcoord_pixel=(Xcoord-membrane_tl[1])/resolution_x+120, #120 added because there seem to be a consistent offset
         Ycoord_pixel=(Ycoord-membrane_tl[2])/resolution_y)

microsample_centroid_pixel <- input_data %>% 
  summarise(mean_x=mean(Xcoord_pixel),
            mean_y=mean(Ycoord_pixel))

input_data %>% 
  ggplot(aes(x=Xcoord_pixel,y=-Ycoord_pixel,color=SampleType)) +
    geom_point()

#Visual confirmation

# Crop top-left coordinates
crop_ref_x=round(microsample_centroid_pixel$mean_x-(width/2),0)
crop_ref_y=round(microsample_centroid_pixel$mean_y-(height/2),0)

# Adjust top-left if it's over the edge
if(crop_ref_x<20){crop_ref_x=20}
if(crop_ref_y<20){crop_ref_y=20}


# Load input figure
input_image <- image_read(overview_image)

#Crop input figure and save
cropped_image <- image_crop(input_image, geometry = sprintf("%dx%d+%d+%d", width, height, crop_ref_x, crop_ref_y))
image_write(cropped_image, str_c("output/",cryosection,".jpg"))

# Calculate microsample coordinates in cropped image
input_data <- input_data %>% 
  mutate(Xcoord_pixel_crop=round(Xcoord_pixel-crop_ref_x,0), #120 added because there seem to be a consistent offset
         Ycoord_pixel_crop=round(Ycoord_pixel-crop_ref_y,0))


# Print microsamples in cropped image
cropped_image2 <- image_draw(cropped_image)
for (i in c(1:nrow(input_data))) {
  cropped_image2 <- image_draw(cropped_image2)
  points(input_data$Xcoord_pixel_crop[i], input_data$Ycoord_pixel_crop[i], col = "red", pch = 16, cex = 2)  # Change `cex` for dot size
  dev.off()  # Finish drawing
}
image_write(cropped_image2, str_c("output/",cryosection,"_marked.jpg"))

# Output cryosection csv
input_data %>% 
  select(microsample,cryosection,SampleType,Xcoord,Ycoord,size,shape,Xcoord_pixel_crop,Ycoord_pixel_crop) %>% 
  write_csv(str_c("output/",cryosection,".csv"))
