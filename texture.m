

function image=texture(input)

input=imresize(input,[65 65]);



%%%%%%%%%%%dividing image into 16x16 blocks%%%%%%%%%%%
[rows ,columns ,numberOfColorBands] = size(input);
blockSizeR = 5; % Rows in block.
blockSizeC = 5; % Columns in block.       
wholeBlockRows = rows / blockSizeR;
blockVectorR = [blockSizeR * ones(1, wholeBlockRows)];
wholeBlockCols = columns / blockSizeC;
blockVectorC = [blockSizeC * ones(1, wholeBlockCols)];

if numberOfColorBands > 1
	% It's a color image.
	ca = mat2cell(input, blockVectorR, blockVectorC, numberOfColorBands);
    M=mat2cell(input, blockVectorR, blockVectorC, numberOfColorBands);
else
	ca = mat2cell(input, blockVectorR, blockVectorC);
    M=mat2cell(input, blockVectorR, blockVectorC);
end
% Now display all the blocks.


plotIndex = 1;
numPlotsR = size(ca, 1);
numPlotsC = size(ca, 2);




for r = 1 : numPlotsR
	for c = 1 : numPlotsC
        figure(2),
		subplot(numPlotsR, numPlotsC, plotIndex);
		rgbBlock = ca{r,c};
		imshow(rgbBlock);
       
        
        
         plotIndex = plotIndex + 1;
     
    end
 
 

end



k=1;
for i=1:13
 for j=1:13
    for m=1:5
        for n=1:5
          
        Y=ca{i,j}(m,n)/(5*5);
        delx{i,j}(m,n)=(ca{i,j}(m,n)-Y);
        
      
        end
    end
    
      figure(3),
		subplot(13,13,k);
        imshow(delx{i,j})
        
       
 
      
    k=k+1; 
     
 end    
end






 image=cell2mat(delx);
 figure,imshow(image)

 
 ima1=255-image;
 figure,imshow(ima1,[])
   
end


