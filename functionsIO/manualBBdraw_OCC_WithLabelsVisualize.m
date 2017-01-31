function modifiedImage=manualBBdraw_OCC_WithLabelsVisualize(img,bb,bbOCC,myColor,myColorOCC,lWidth,myText1,myText2,myFigNumber)

% MANUALBBDRAW_OCC_WITHLABELSVISUALIZE.m produces graphical output for DSKCF tracker
%
%  MANUALBBDRAW_OCC_WITHLABELSVISUALIZE overlays on top of the current
%  frame the tracker bounding box and the one corresponding to the
%  occluding object (if they exist) with different color and linewidth.
%  Also a text is displayed on the top left corner of the image
%
%   INPUT:
%  -img source image
%  -bb target bounding box bounding box in the format [topLeftX, topLeftY,
%   bottomRightX, bottomRightY] read as [columnIndexTopLeft, rowIndexTopLeft,
%   columnIndexBottomRight, rowIndexBottomRight]
%
%  -bbOCC occluding object bounding box bounding box in the format
%  [topLeftX, topLeftY,bottomRightX, bottomRightY] read as
%  [columnIndexTopLeft,rowIndexTopLeft, columnIndexBottomRight,
%  rowIndexBottomRight]
%
%  -myColor and myColorCC, selected color for the bounding boxes
%  -lWidth selected line width for the two bounding boxes
%  -myText1,myText2 text related to the tracker and occluding object
%  displayed (with the selected color) on the top left of the image
%
%  See also ZBUFFER_CDATA, MANUALBBDRAW_OCC_WITHLABELS, MANUALBBDRAW_OCC
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk

imageSize=[size(img,2),size(img,1)];

figure(myFigNumber)
%myFig=figure('visible','off');
%set(myFig,'resize','off');
hold on;
imshow(img);
if(isempty(myText1)==false)
    text(25,25,myText1,'color',myColor,'fontsize',20,'fontweight','bold');
end
if(isempty(myText2)==false)
    text(25,60,myText2,'color',myColorOCC,'fontsize',20,'fontweight','bold');
end


if(isempty(bb)==false & isnan(bb)==false)
    if(bb(1)>imageSize(1) | bb(2)>imageSize(2))
        bb=[];
    else
        bb(bb(1:2)<0)=1;
        if(bb(1)+bb(3)>imageSize(1))
            bb(3)=imageSize(1)-bb(1);
        end
        
        if( bb(2)+bb(4)>imageSize(2))
            bb(4)=imageSize(2)-bb(2);
        end
        
        rectangle('Position', bb,'LineWidth',lWidth,'edgecolor',myColor);
    end
end



if(isempty(bbOCC)==false & isnan(bbOCC)==false)
    
    if(bbOCC(1)>imageSize(1) | bbOCC(2)>imageSize(2))
        bbOCC=[];
    else
        bb(bbOCC(1:2)<0)=1;
        if(bbOCC(1)+bbOCC(3)>imageSize(1))
            bbOCC(3)=imageSize(1)-bbOCC(1)+1;
        end
        
        if( bbOCC(2)+bbOCC(4)>imageSize(2))
            bbOCC(4)=imageSize(2)-bbOCC(2)+1;
        end
        rectangle('Position', bbOCC,'LineWidth',lWidth,'edgecolor',myColorOCC,'lineStyle','--');
    end
end

F = im2frame(zbuffer_cdata(gcf));
modifiedImage=imresize(F.cdata,[size(img,1),size(img,2)]);

%close (myFig)
%pause()