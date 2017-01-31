function modifiedImage=manualBBdraw_MULTIPLE(img,bbVector,myColor,lWidth,myText1,myTextColor)

% MANUALBBDRAW_MULTIPLE.m function that produces graphical output for DSKCF tracker
%
%  MANUALBBDRAW_MULTIPLE overlays on top of the current frame a
%  vector of bounding boxes (if they exist) with different color and
%  linewidth. Also a text is displayed on the top left corner of the image
%
%   INPUT:
%  -img source image
%  -bbVector bounding box vector in the format [topLeftX, topLeftY,
%   bottomRightX, bottomRightY] read as [columnIndexTopLeft, rowIndexTopLeft,
%   columnIndexBottomRight, rowIndexBottomRight]
%
%
%  -myColor selected color for the bounding boxes
%  -lWidth selected line width for the bounding boxes
%  -myText1 list containing text corresponding to the bounding boxes
%  displayed (with the selected color) on the top left of the image
%
%   OUTPUT
%  -modifiedImage captured from matlab fig, ready to be saved
%
%  See also MANUALBBDRAW_OCC, MANUALBBDRAW_OCC_WITHLABELSVISUALIZE, ZBUFFER_CDATA
%
%  University of Bristol
%  Massimo Camplani and Sion Hannuna
%
%  massimo.camplani@bristol.ac.uk
%  hannuna@compsci.bristol.ac.uk

imageSize=[size(img,2),size(img,1)];


myFig=figure('visible','off');
set(myFig,'resize','off');
hold on;
imshow(img);
if(isempty(myText1)==false)
    if(iscell(myText1)==0)
        text(25,25,myText1,'color',myColor(1),'fontsize',20,'fontweight','bold');
    else
        for jj=1:length(myText1)
            text(25,25+35*(jj-1),myText1(jj),'color',myTextColor(jj),'fontsize',20,'fontweight','bold');
        end
    end
end

for kk=1:size(bbVector,1)
bb=bbVector(kk,:);
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
        
        rectangle('Position', bb,'LineWidth',lWidth,'edgecolor',myColor(kk));
    end
end

end


F = im2frame(zbuffer_cdata(gcf));

modifiedImage=imresize(F.cdata,[size(img,1),size(img,2)]);

close (myFig)
