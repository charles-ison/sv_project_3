# Project 3
Scientific Visualizations  
Charles Ison  
11/13/2022

GitHub Link (for README markdown formatting): https://github.com/charles-ison/sv_project_3

### Running The Program
* Project was compiled and tested using VisualStudio on the machines in the Kelly computer lab  
* To change the photo that is used for each answer, update the file path in the learnply.cpp keyboard() (images are stored in the data/images/... folder)
* To run part 1.a.1 use the 'a' button (displays chosen image as input to the IBFV algorithm)
* To run part 1.a.2 use the 'b' button (displays chosen image as input to the IBFV algorithm and uses the vector field)
* To run part 1.c. use the 'c' button (extracts edges from the chosen image using the Sobel filter)
* To run part 1.d. use the 'd' button (extracts edges from the chosen image using the Sobel filter and uses the vector field)
* To run part 2.a. use the 'e' button (extracts all singularities from the vector field, the key to match the color and type is below)
* To run part 2.b. there are two ways to select a point, first either "ctrl" + "click" on a quad's vertex and a forward and backward streamline will be generated or "shift" + "click" on a quad and a forward and backward streamline will be generated from vertex V2
* To run part 3.a. use the 'f' button (extracts all singularities from the vector field on top of the IBFV vizualization of the same field)
* To run part 3.b. use the 'g' button (extracts all singularities from the vector field on top a streamline vizualization)
* To run part 3.c. use the 'h' button (extracts all singularities and separatrices from the vector field on top of a streamline visualization)

### Singularity Color Key
* Red is a source
* Blue is a sink
* Green is a saddle
* Light Blue is a center
* Yellow is a focus

### Question Answers

1.a.) When we use the image as the initial noise into the IBFV algorithm, the algorithm seems to smear to image in the directions the vector field is going. It creates a almost dreamy effect where the image starts to feel surreal. The impact can vary a lot depending on both the vector field and the image, but for the vector field visualized below (V10), I liked how it caused the image to almost spiral inward framing the mountains in each images.

* Avalanche Lake - Glacier National Park
![image](mountains1_part1a.png)
* Bends the trees that frame the mountain image, almost feels like a Salvador Dali painting to me.
![image](mountains1_part1b.png)

* Driving - Grand Teton National Park
![image](mountains2_part1a.png)
* Makes it seem a bit like you are driving along the road towards the mountains at warp speed.
![image](mountains2_part1b.png)

* Cracker Lake - Glacier National Park
![image](mountains3_part1a.png)
* There are already patterns and lines in the rock face that blend in with the vector field. It is not immediately clear which parts of the rock face are natural or caused by the vector field.
![image](mountains3_part1b.png)

1.b.) I personally like the V10 vector field the most for each of the three mountain lakes. First, I felt like the V10 vector fields creates a more dramatic change to the image than many of the other vector fields and I just enjoy the way while spiraling inward, it creates both a central focus area and a framing around the outside of the mountains.

1.c.)

* Avalanche Lake - Glacier National Park
* I really enjoy this visualization, I think the contrast between the trees, the lake, the waterfalls and the cliff become more noticeable than in the normal image.
![image](mountains1_part1c.png)

* Driving - Grand Teton National Park
* In this one the road remains clear, but it becomes harder to tell what the mountains actually are. This could partically be do to the poor initial image quality after compressing.
![image](mountains2_part1c.png)

* Cracker Lake - Glacier National Park
* In this one, the lake becomes very clear and looks calm, but the rock face is very busy.
![image](mountains3_part1c.png)

1.d.)

* Avalanche Lake - Glacier National Park
* Compared to part 1.b., I think I enjoy the edge field visualizaiton more due to the contrast it makes between the main image photos: water, trees, and mountains. Once again, the V10 vector field was my favorite, but even more so on the edge field because it creates interesting framing with the bending trees versus the lake in the center.
![image](mountains1_part1d.png)

* Driving - Grand Teton National Park
* Compared to part 1.b., I also like this image more. One reason for this is I feel the overall image quality is a bit low which becomes more noticeable in part 1.b., but here because we are focusing on the edges it seems more clear. Also the vector field V10 makes it seems like the mountains are spiraling inward while you are driving towards them, which is just a unique visual effect.
![image](mountains2_part1d.png)

* Cracker Lake - Glacier National Park
* Compared to part 1.b., this might be the only one I liked less. Since the rock face is such a big part of the photo, it becomes harder to tell what is actually happening in 1.d. and I think the colors are unique and have a big impact in 1.b. That being said, the V10 vector field was once again my favorite, due to its dreamy and surreal impact on the photo.
![image](mountains3_part1d.png)

2.a.)

* V1
* See "Singularity Color Key" above for classification
![image](v1_2a.png)

* V3
* See "Singularity Color Key" above for classification
![image](v3_2a.png)

* V4
* See "Singularity Color Key" above for classification
![image](v4_2a.png)

* V5
* See "Singularity Color Key" above for classification
![image](v5_2a.png)

* V6
* See "Singularity Color Key" above for classification
![image](v6_2a.png)

* V8
* See "Singularity Color Key" above for classification
![image](v8_2a.png)

* V9
* See "Singularity Color Key" above for classification
![image](v9_2a.png)

* V10
* See "Singularity Color Key" above for classification
![image](v10_2a.png)

2.b.)

* V1
* See "Running The Program" for directions on how to click and generate a single forward and backward streamline.
![image](v1_2b.png)

* V3
* See "Running The Program" for directions on how to click and generate a single forward and backward streamline.
![image](v3_2b.png)

* V4
* See "Running The Program" for directions on how to click and generate a single forward and backward streamline.
![image](v4_2b.png)

* V5
* See "Running The Program" for directions on how to click and generate a single forward and backward streamline.
![image](v5_2b.png)

* V6
* See "Running The Program" for directions on how to click and generate a single forward and backward streamline.
![image](v6_2b.png)

* V8
* See "Running The Program" for directions on how to click and generate a single forward and backward streamline.
![image](v8_2b.png)

* V9
* See "Running The Program" for directions on how to click and generate a single forward and backward streamline.
![image](v9_2b.png)

* V10
* See "Running The Program" for directions on how to click and generate a single forward and backward streamline.
![image](v10_2b.png)

3.a.)

* V1
* The extracted singularity matches the IBFV visualization as it is a focus on the IBFV visualization shows the vector field spiraling out from the singularity.
![image](v1_3a.png)

* V3
* The extracted singularity matches the IBFV visualization as it is a saddle and the IBFV visualization shows the vector field flowing around the singularity.
![image](v3_3a.png)

* V4
* The extracted singularity matches the IBFV visualization as it is a saddle and the IBFV visualization shows the vector field flowing around the singularity.
![image](v4_3a.png)

* V5
* The extracted singularities mostly match the IBFV visualization as there are two focuses which the IBFV visualization shows the vector field spiraling out from both. There is also a single sink, which might be right be right or wrong, but I can not clearly enough see from the IBFV visualization to tell.
![image](v5_3a.png)

* V6
* The extracted singularities match the IBFV visualization as there is a green saddle the IBFV visualization shows the vector field flowing around and a yellow focus the IBFV visualization shows the vector field spiraling out from.
![image](v6_3a.png)

* V8
* The extracted singularities match the IBFV visualization as there is a green saddle the IBFV visualization shows the vector field flowing around and a yellow focus the IBFV visualization shows the vector field spiraling out from.
![image](v8_3a.png)

* V9
* There are seven extracted singularities and it is hard to tell if they all match the IBFV visualization, but I suspect some might be wrong. I feel I can definitely see two saddles and two focuses in the IBFV visualization, which match four of our extracted singularities, but it is hard to tell for the other singularities without a higher resolution image and the ability to zoom in. One reason there could be incorrectly extracted singularities is we rely on an epsilon threshold for extracting the singularities which might need further tuning.
![image](v9_3a.png)

* V10
* There are four extracted singularities, two sources and two focuses. Much like V9, I have a hard time telling 100% if the extracted singularities match the IBFV visualization here. But, the upper source and upper focus do look visually correct to me. Where the bottom source and bottom focus are clustered, I do think there is a singularity there, but I cannot tell the type or the amount from the visualization. One reason there could be incorrectly extracted singularities is we rely on an epsilon threshold for extracting the singularities which might need further tuning.
![image](v10_3a.png)

3.b.)

* V1
* One strength of the streamline here is I feel it more clearly makes the flow of the vector field around singularity more clear. A weakness of the streamlines here versus IBFV is since we are manually generating them iteratively, it is harder to see what is happening in all areas of the field as they streamlines tend to cluster. One other weakness is the streamlines tend to get more faint toward singularities, I am not clear at this moment why this occurs.
![image](v1_3b.png)

* V3
* Same strengths and weaknesses mentioned earlier, but this time almost half the screen does not contain any information about the vector fields topology due to the previously mentioned weakness.
![image](v3_3b.png)

* V4
* Same strengths and weaknesses mentioned earlier, but this time more than half the screen does not contain any information about the vector fields topology due to the previously mentioned weakness.
![image](v4_3b.png)

* V5
* Same strengths and weaknesses mentioned earlier, but this time almost half the screen does not contain any information about the vector fields topology due to the previously mentioned weakness.
![image](v5_3b.png)

* V6
* Same strengths and weaknesses mentioned earlier, but this time almost half the screen does not contain any information about the vector fields topology due to the previously mentioned weakness.
![image](v6_3b.png)

* V8
* Same strengths and weaknesses mentioned earlier, but this time most of the screen is conveying information about the topology, so the previous weakness is not as bad.
![image](v8_3b.png)

* V9
* Same strengths and weaknesses mentioned earlier, but this time more than half the screen does not contain any information about the vector fields topology due to the previously mentioned weakness. Also, there are very few streamlines near the singularities, so there is very little information given to know their existance or classification just based on the streamlines.
![image](v9_3b.png)

* V10
* Same strengths and weaknesses mentioned earlier, but this time almost half the screen does not contain any information about the vector fields topology due to the previously mentioned weakness. Unlike V9, there are more streamlines around the singularities, but still not as much information is conveyed as the IBFV algorithm. 
![image](v10_3b.png)

3.c.)

* V1
* No separatrices are found, so the visualization is essentionall the same as part 3.b.
![image](v1_3c.png)

* V3
* Two eparatrices are found and help give further evidence that the singularitiy is a saddle. The visualization provides more overall info about the singularities than part 3.b.
![image](v3_3c.png)

* V4
* Two eparatrices are found and help give further evidence that the singularitiy is a saddle. The visualization provides more overall info about the singularities than part 3.b.
![image](v4_3c.png)

* V5
* No separatrices are found, so the visualization is essentionall the same as part 3.b.
![image](v5_3c.png)

* V6
* Two separatrices are found and gives further evidence to suppor the saddle and focus singularity classifications. The visualization provides more overall info about the singularities than part 3.b.
![image](v6_3c.png)

* V8
* Two separatrices are found and gives further evidence to suppor the saddle and focus singularity classifications. The visualization provides more overall info about the singularities than part 3.b.
![image](v8_3c.png)

* V9
* Many separatrices are found and give much more information about the singurities than what we were able to visualize just using streamlines. Many of the seven singularities were borderline impossible to confirm their existance and classify using just the streamlines, but now all but one out of the seven seem to be correctly classified. The lower focus still does not have enough information to tell 100%.
![image](v9_3c.png)

* V10
* No separatrices are found, so the visualization is essentionall the same as part 3.b.
![image](v10_3c.png)
