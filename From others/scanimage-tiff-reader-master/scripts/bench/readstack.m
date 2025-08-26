function v=readstack(fname)

    info=imfinfo(fname);
    im=imread(fname,1);
    v=zeros([size(im) numel(info)],'like',im);
    v(:,:,1)=im;
    for i=2:numel(info)
        v(:,:,i)=im;
    end
%{
    t=Tiff(fname,'r');
    sz=[t.getTag('ImageWidth') t.getTag('ImageLength') nifd];
    t.setDirectory(1);
    im=t.read;
    v=zeros(sz,'like',im);
    v(:,:,1)=im;
    for i=2:nifd
        t.setDirectory(i)
        v(:,:,i)=t.read;
    end
    
    function n=nifd
        t.setDirectory(1);
        try
            while 1,
                t.nextDirectory;
            end
        catch
            
        end
        n=t.currentDirectory;
    end
%}
end