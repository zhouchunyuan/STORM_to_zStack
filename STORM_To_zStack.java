import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.util.ArrayUtil;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.measure.*;
import ij.text.*;
import ij.io.*;
import java.io.PrintWriter;
import java.io.*;
import java.io.FileReader;
import java.util.ArrayList;
import ij.plugin.filter.Filler;
import java.awt.event.*;
import java.util.Arrays;
import java.util.List;
import javax.swing.JOptionPane;

/////////////////////////////////////////////////////////
// 3D-STORM to Zstack tool
/////////////////////////////////////////////////////////
public class STORM_To_zStack implements PlugIn {

    int totalFrames;//
    boolean keepZreject = false;
    boolean showTable = false;
    int windowSize = 512;

    public String[] setParameters(){
        String strfiles="";

        OpenDialog openDlg = new OpenDialog("Choose a STORM particle tables","","*.txt");
        //filename = openDlg.getPath();
        ///////
        if(openDlg.getDirectory()==null)return null;
        File[] fileList;
        File dir = new File(openDlg.getDirectory());

        fileList= dir.listFiles(new FilenameFilter() { 
                public boolean accept(File dir, String filename)
                { return filename.endsWith(".txt"); }
                });

        boolean[] fileSelected=new boolean[fileList.length];
        ///////
        GenericDialog gd = new GenericDialog("Input parameters");
        String [] tmp= new String[fileList.length];
        for(int i=0;i<fileList.length;i++){
            tmp[i]=fileList[i].getName();
            if(fileList[i].getPath().equals(openDlg.getPath()))fileSelected[i]=true;
        }
        gd.addMessage("STORM files:");
        gd.addCheckboxGroup(fileList.length,1,tmp,fileSelected);
        //gd.addStringField("Title: ", title);
        gd.addMessage("Options:");
        gd.addCheckbox("Show \"z rejected dots\"",keepZreject);
        gd.addCheckbox("Show STORM table",showTable);

        gd.showDialog();
        if (gd.wasCanceled()) return null;
        //title = gd.getNextString();
        for(int i=0;i<fileList.length;i++){
            fileSelected[i]=gd.getNextBoolean();
            if(fileSelected[i])strfiles+=fileList[i].getPath()+",";
        }
        keepZreject = gd.getNextBoolean();
        showTable = gd.getNextBoolean();
        return strfiles.split(",");
    }

    // program main()
    // init of 2 tables
    public void run(String arg){
        if (IJ.versionLessThan("1.48v")) return;

        String [] tableFiles;
        do{
            tableFiles = setParameters();
            if(tableFiles==null)return;
            if(tableFiles.length!=1)
                IJ.error("You must select only one STORM table !");
        }while(tableFiles.length!=1);
        ImageTable itbl = new ImageTable();
        itbl.filename = tableFiles[0];
        itbl.run();
        //if(neighborCorrection){
        //}
    }

    ///////////////////////////////////////////////////////
    // This is the main class.
    //
    // It contains the ResultsTable for all dots.
    // All dots are also copied into an NxM dimention matrix "chainArray".
    // Each element in the chainArray contains the dots inside this area.
    //
    // The identified dot-groups are stored in ArrayList<groupInfo>
    // The image window can handle mouse events
    ////////////////////////////////////////////////////////
    public class ImageTable {

        int fileID;
        String filename;
        ImagePlus imp;
        ImagePlus []impStk;
        Overlay overlay = new Overlay();//image overlay is used to compare Master/Slave images
        ResultsTable tbl;
        ResultsTable corr;//correction of xy
        List<String> channelNamesArray = new ArrayList<String>();

        // will contain displacement in xyz
        double dispX = 0.0;
        double dispY = 0.0;
        double dispZ = 0.0;

        //z range
        double zMin,zMax;

        //maxPhoton
        double maxPhotons=0;

        double zoom=0.1;// means pixels/nm
        double x_L; 
        double x_R;
        double y_U;
        double y_D;

        // below is a NxM matrix area map 
        // each area has many dots connected by 'DotChain'
        int N=8;
        int M=8;
        //set full size to be 120% of 40960 nm
        double fullWidth = 4096*10;// nm
        double fullHeight = 4096*10;// nm image.setWindow(new Window(image, image.getCanvas())); 
        double matrixSizeX = fullWidth/N;
        double matrixSizeY = fullHeight/M;

        //one single dot used by ArrayList<dot> Dots
        private class DotChain {
            int ch;
            double x;    
            double y;
            double z;
            double photons;
            double accuracy;
            int frame;
            DotChain next=null;
            DotChain(double a, double b, double c, int n){
                this(a,b,c);
                frame = n;
            }
            DotChain(double a, double b, double c){
                x=a;y=b;z=c;
            }
            DotChain(double a, double b, double c, int n,double d, double e){
                this(a,b,c,n);
                photons=d;
                accuracy=e;
            }
            DotChain(double a, double b, double c, int n,double d, double e, int channel){
                ch=channel;
                x=a;y=b;z=c;
                frame = n;
                photons=d;
                accuracy=e;
            }


            void attach(DotChain dot){
                if (this.next==null)
                    this.next=dot;
                else
                    this.next.attach(dot);
            }
        }//Class DotChain

        DotChain[][] chainArray=new DotChain[N][M];

        private void clearDotChain(){
            for(int i=0;i<N;i++)
                for(int j=0;j<M;j++)
                    chainArray[i][j]=null;
        }

        public void ImageTable() {}

        public void run(){
            //setParameters();
            //OpenDialog openDlg = new OpenDialog("Choose an N-STORM particle table","","*.txt");
            //filename = openDlg.getPath();
            readFromFile(filename);
            fillChainArray();

            //getCorrectionTable();
            imp= NewImage.createByteImage("STORM window ["+fileID+"]",windowSize,windowSize,1,NewImage.FILL_BLACK);

            imp.show();

            CustomCanvas cc = new CustomCanvas(imp,this);
            imp.setWindow(new ImageWindow(imp,cc)); 

            zoomInWindow(0.1,fullWidth/2,fullHeight/2);

        }

        //fill the chain array using "tbl"
        private void fillChainArray(){
            // DotChain[][] chainArray=new DotChain[N][M];
            clearDotChain();
            int chName_col = tbl.getColumnIndex("Channel Name");
            int xc_col = tbl.getColumnIndex("Xc");
            int yc_col = tbl.getColumnIndex("Yc");
            int zc_col = tbl.getColumnIndex("Zc");
            int frame_col = tbl.getColumnIndex("Frame");
            int acc_col = tbl.getColumnIndex("Lateral Localization Accuracy");
            int pho_col = tbl.getColumnIndex("Photons");
            int count = tbl.getCounter();

            for(int i=0;i<count;i++){
                int ch=0;
                String chName = tbl.getStringValue(chName_col,i).trim();
                for(int chNum = 1;chNum<channelNamesArray.size();chNum++)
                    if(chName.equals(channelNamesArray.get(chNum).trim()))ch=chNum; 
                double x = tbl.getValueAsDouble(xc_col,i);
                double y = tbl.getValueAsDouble(yc_col,i);
                double z = tbl.getValueAsDouble(zc_col,i);
                double pho = tbl.getValueAsDouble(pho_col,i);//photons
                if(pho>maxPhotons)maxPhotons=pho;//find max photon number
                double acc = tbl.getValueAsDouble(acc_col,i);//lateraal accuracy
                int frame =(int) tbl.getValueAsDouble(frame_col,i);
                int n = (int)(x/matrixSizeX);
                int m = (int)(y/matrixSizeY);
                if(n==N)n--;
                if(m==M)m--;
                if (x>fullWidth || y>fullHeight){
                    IJ.log("("+i+") x="+x+", y="+y+" skipped");
                    continue;
                }
                if(chainArray[n][m]==null)
                    chainArray[n][m]=new DotChain(x,y,z,frame,pho,acc,ch);
                else
                    chainArray[n][m].attach(new DotChain(x,y,z,frame,pho,acc,ch));
                IJ.showProgress(i,count);
            }

        }//function fillChainArray()



        //create correction table using current particle table
        // be careful: there may be missing frames!!
        private void getCorrectionTable(){
            int x_col = tbl.getColumnIndex("X");
            int y_col = tbl.getColumnIndex("Y");
            int z_col = tbl.getColumnIndex("Z");
            int xc_col = tbl.getColumnIndex("Xc");
            int yc_col = tbl.getColumnIndex("Yc");
            int zc_col = tbl.getColumnIndex("Zc");
            int frame_col = tbl.getColumnIndex("Frame");
            int count = tbl.getCounter();
            int frame_old = -1;// old frame index, judge of frame change
            corr = new ResultsTable();

            for(int i=0;i<count;i++){
                double dx = -(tbl.getValueAsDouble(xc_col,i)- tbl.getValueAsDouble(x_col,i));
                double dy = -(tbl.getValueAsDouble(yc_col,i)- tbl.getValueAsDouble(y_col,i));
                double dz = -(tbl.getValueAsDouble(zc_col,i)- tbl.getValueAsDouble(z_col,i));
                int frame =(int) tbl.getValueAsDouble(frame_col,i);
                if(frame_old != frame){
                    corr.incrementCounter();
                    corr.addValue("frame",frame);
                    corr.addValue("dx",dx);
                    corr.addValue("dy",dy);
                    corr.addValue("dz",dz);
                    frame_old = frame;
                }
            }
            //corr.show("correction table"+fileID);

        }//function getCorrectionTable()

        //read particle table and copy into ResultsTable tbl
        private void readFromFile(String path){
            tbl = new ResultsTable();
            try{
                int idx,arrayLen;
                String lineString;
                FileReader fr = new FileReader(path);
                BufferedReader br = new BufferedReader(fr);
                String [] labelElements,elements;
                int counter = 1;

                // read the labels
                labelElements = br.readLine().split("\t");
                arrayLen = labelElements.length;

                // Checking an String in String Array by converting Array To ArrayList
                List<String> labelList = Arrays.asList(labelElements);
                //int col_chName = labelList.indexOf("Channel Name");
                int col_Xc = labelList.indexOf("Xc");
                int col_Yc = labelList.indexOf("Yc");
                int col_Zc = labelList.indexOf("Zc");// used to find z range
                int col_Frame = labelList.indexOf("Frame");

                //fill out the particle table
                while ((lineString = br.readLine()) != null)
                {
                    IJ.showStatus("Loading N-STORM particle data");
                    counter++;
                    lineString.trim();
                    elements = lineString.split("\t");
                    if(elements[0].trim().equals("Z Rejected") && !keepZreject)
                    {continue;} 
                    if(!channelNamesArray.contains(elements[0]))channelNamesArray.add(elements[0]);
                    tbl.incrementCounter();
                    tbl.addValue(labelElements[0],elements[0]);
                    //ctable.addLabel(labelElements[0],elements[0]);
                    // This will add a string as label, which is not counted as values
                    for( idx = 1;idx<arrayLen;idx++)
                    {
                        tbl.addValue(labelElements[idx],Double.parseDouble(elements[idx]));
                    }
                    //re-dim full size
                    if(fullWidth<Double.parseDouble(elements[col_Xc]))
                        fullWidth=Double.parseDouble(elements[col_Xc]);
                    if(fullHeight<Double.parseDouble(elements[col_Yc]))
                        fullHeight=Double.parseDouble(elements[col_Yc]);

                    //only last change is useful
                    totalFrames =(int)Double.parseDouble(elements[col_Frame] );
                }
                IJ.showStatus("N-STORM data loaded. Creating image window...");
                //re-dim matrice size
                matrixSizeX = fullWidth/N;
                matrixSizeY = fullHeight/M;
                fr.close();
                //find the zMin and zMax
                float [] zArray = tbl.getColumn(col_Zc);
                ArrayUtil aUtil = new ArrayUtil(zArray);
                zMin = aUtil.getMinimum();
                zMax = aUtil.getMaximum();
                //IJ.log(""+zMin + " ~ " + zMax);
                /*for(int i=0;i<channelNamesArray.size();i++)
                  IJ.log(channelNamesArray.get(i));*/

                if(showTable)tbl.show("N-STORM particle data");

                return ;
            } catch (FileNotFoundException e) {
                IJ.error("File not found exception" + e);
                return ;
            } catch (IOException e) {
                IJ.error("IOException exception" + e);
                return ;
            } catch (NumberFormatException e) {
                IJ.error("Number format exception" + e);
                return ;
            }
        }//function readFromFile()

        //display zoomed particles
        private void zoomInWindow(double zoom, double x0, double y0){
            this.zoom = zoom;
            int imgWidth = imp.getWidth();
            int imgHeight= imp.getHeight();
            x_L= x0-imgWidth/2/zoom;
            x_R= imgWidth/2/zoom+x0;
            y_U= y0-imgHeight/2/zoom;
            y_D= imgHeight/2/zoom+y0;
            //IJ.log("x_L="+x_L+"x_R="+x_R+"y_U="+y_U+"y_D="+y_D);
            ImageProcessor ip = imp.getProcessor(); 

            //Filler filler = new Filler();
            //filler.clear(ip); //this in some cases will fill into white color.
            ip.setColor(0);
            ip.fill();

            ip.setLineWidth(2);//set the drawDot style
            ip.setColor(255);
            //frame should shift by dispX,dispY
            int n1 = (int)Math.floor((x_L-dispX)/matrixSizeX);
            if (n1<0) n1=0;
            int n2 = (int)Math.ceil((x_R-dispX)/matrixSizeX);
            if (n2>N-1)n2=N-1;
            int m1 = (int)Math.floor((y_U-dispY)/matrixSizeY);
            if (m1<0) m1=0;
            int m2 = (int)Math.ceil((y_D-dispY)/matrixSizeY);
            if (m2>M-1) m2 =M-1;
            for(int n = n1; n<=n2; n++)
                for(int m = m1; m<=m2; m++){
                    if(chainArray[n][m]!=null){
                        DotChain dot=chainArray[n][m];
                        while(true){
                            int x = (int)((dot.x+dispX-x_L)*zoom);//dispX added
                            int y = (int)((dot.y+dispY-y_U)*zoom);//dispY added
                            ip.drawDot(x,y);
                            //ip.putPixel(x,y,255);
                            //IJ.log("("+dot.x+","+dot.y+")=>("+x+","+y+")");
                            if(dot.next==null)
                                break;
                            else
                                dot=dot.next;
                        }
                    }
                }

            IJ.run(imp, "Set Scale...", "distance="+zoom+" known=1 pixel=1 unit=nm");

            imp.updateAndDraw();

        }//function zoomInWindow


        public double findLocalMaxPhotons(){
            double max=0;
            int n1 = (int)Math.floor((x_L-dispX)/matrixSizeX);
            if (n1<0) n1=0;
            int n2 = (int)Math.ceil((x_R-dispX)/matrixSizeX);
            if (n2>N-1)n2=N-1;
            int m1 = (int)Math.floor((y_U-dispY)/matrixSizeY);
            if (m1<0) m1=0;
            int m2 = (int)Math.ceil((y_D-dispY)/matrixSizeY);
            if (m2>M-1) m2 =M-1;
            for(int n = n1; n<=n2; n++)
                for(int m = m1; m<=m2; m++){
                    if(chainArray[n][m]!=null){
                        DotChain dot=chainArray[n][m];
                        while(true){
                            int x = (int)((dot.x+dispX-x_L)*zoom);//dispX added
                            int y = (int)((dot.y+dispY-y_U)*zoom);//dispY added
                            if(x>=0 && x<windowSize)
                                if(y>=0 && y<windowSize)
                                    if(dot.photons>max)max=dot.photons;
                            if(dot.next==null)
                                break;
                            else
                                dot=dot.next;
                        }
                    }
                }
            return max;

        }//end of findLocalMaxPhotons()

        public void createStack(){
            double zstep = 1/zoom;
            double dotFactor = 4.0;//zoom factor for Gaussian dot
            // double psize = 1.0;
            boolean useLocalMaxPhotons=true;
            boolean showMerge=false;
            boolean showInfo = true;
            double localMaxPhotons;
            int channels = channelNamesArray.size();
            String [] radioItems = {"3D Gaussian","2D Gaussian","3D Ball"};
            int radioChoice = 0;
            int nProcessed = 0;

            GenericDialog gd = new GenericDialog("Create Stack From 3D STORM data");
            gd.addMessage("Options:");
            gd.addMessage("Will create "+channels+" channels");
            gd.addCheckbox("Brightness control by local max photon number",useLocalMaxPhotons);

            gd.addRadioButtonGroup("Choose a rendering method:",radioItems,radioItems.length,1,"3D Gaussian");
            gd.addMessage("(notes:the Z step is recommended to be equal to XY (="+zstep+")");
            gd.addNumericField("Z Step:",zstep,2,5,"nm");
            gd.addSlider("Dot size factor:", 1.0, 10.0, dotFactor);
            gd.addCheckbox("Show merge image",showMerge);
            gd.addCheckbox("Show info in log",showInfo);

            gd.showDialog();
            if (gd.wasCanceled()) return ;
            //title = gd.getNextString();
            useLocalMaxPhotons = gd.getNextBoolean();
            zstep = gd.getNextNumber();
            dotFactor = gd.getNextNumber();
            String item = gd.getNextRadioButton();
            showMerge = gd.getNextBoolean();
            showInfo = gd.getNextBoolean();
            //return strfiles.split(",");

            for(int i=0;i<radioItems.length;i++){
                if(item.equals(radioItems[i]))radioChoice=i;
            }
            //IJ.log(""+radioChoice);

            int zNumber = (int)((zMax-zMin)/zstep);
            impStk = new ImagePlus[channels];
            for(int i=0;i<channels;i++){
                impStk[i]= NewImage.createByteImage(channelNamesArray.get(i)+" STORM Stack ",
                        windowSize,windowSize,zNumber,NewImage.FILL_BLACK);
                Calibration cal = new Calibration();
                cal.pixelHeight = 1/zoom;
                cal.pixelWidth = 1/zoom;
                cal.pixelDepth = zstep;
                cal.setUnit("nm");
                impStk[i].setCalibration(cal);
                //HyperStackConverter hypCvtr = new HyperStackConverter();
                //impStk = hypCvtr.toHyperStack(impStk, 0, zNumber, 0);
                //impStk[i].getProcessor().setColor(Color.white); 
                impStk[i].show();
            }
            //IJ.log("local:"+findLocalMaxPhotons()+" global:"+maxPhotons); 
            //brightness according to local/global maxphotons
            if(useLocalMaxPhotons) 
                localMaxPhotons = findLocalMaxPhotons();
            else
                localMaxPhotons = maxPhotons;

            int n1 = (int)Math.floor((x_L-dispX)/matrixSizeX);
            if (n1<0) n1=0;
            int n2 = (int)Math.ceil((x_R-dispX)/matrixSizeX);
            if (n2>N-1)n2=N-1;
            int m1 = (int)Math.floor((y_U-dispY)/matrixSizeY);
            if (m1<0) m1=0;
            int m2 = (int)Math.ceil((y_D-dispY)/matrixSizeY);
            if (m2>M-1) m2 =M-1;
            for(int n = n1; n<=n2; n++)
                for(int m = m1; m<=m2; m++){
                    if(chainArray[n][m]!=null){
                        DotChain dot=chainArray[n][m];
                        while(true){

                            int x = (int)((dot.x+dispX-x_L)*zoom);//dispX added
                            int y = (int)((dot.y+dispY-y_U)*zoom);//dispY added
                            if(x>=0 && x<windowSize && y>=0 && y<windowSize){
                                int zIdx = (int)((dot.z-zMin)/zstep+0.5);
                                int amp = (int)(255*dot.photons/localMaxPhotons/dotFactor);//amplitude by photons
                                int sigma = (int)(dot.accuracy*zoom*dotFactor);//(zoom is pixel/nm)
                                int ch =dot.ch;
                                nProcessed++;
                                //IJ.log("#"+nProcessed);
                                if(radioChoice == 0)
                                    draw3DGaussian(impStk[ch],x,y,zIdx,sigma,amp);
                                else if(radioChoice == 1)
                                    draw2DGaussian(impStk[ch],x,y,zIdx,sigma,amp);
                                else
                                    draw3DBall(impStk[ch],x,y,zIdx,sigma,amp);

                            }//end of if
                            //ip.putPixel(x,y,255);
                            //IJ.log("("+dot.x+","+dot.y+")=>("+x+","+y+")");
                            if(dot.next==null)
                                break;
                            else
                                dot=dot.next;
                        }
                    }
                }

            if(showMerge){
                try{
                    RGBStackMerge merger = new RGBStackMerge();
                    ImagePlus mergedImp = merger.mergeHyperstacks(impStk,false);
                    mergedImp.show();
                }catch ( Exception ex ){
                    IJ.log("Single channel. Can not create merge!");
                }
            }
            if(showInfo){
                IJ.log(channels+"channel(s) exported into multi-layer images.");
                IJ.log("XY resolution: "+(1/zoom)+" nm");
                IJ.log("Z resolution: "+zstep+" nm");
                IJ.log(nProcessed+" molecules exported");
            }
        }//end of function createStack()

        //draw 2D gaussian
        public void draw2DGaussian(ImagePlus imp,int x0,int y0,int z0,int sigma,int amplitude){
            double fSigma = (double)sigma;
            if(fSigma<0.5)fSigma=0.5;//can not be 0!
            int f = 5;
            ImageProcessor ip = imp.getProcessor(); 
            imp.setSlice(z0);
            for(int y=y0-sigma*f;y<=y0+sigma*f;y++)
                for(int x=x0-sigma*f;x<=x0+sigma*f;x++){
                    double v =Math.exp((-(x-x0)*(x-x0)-(y-y0)*(y-y0))/(2*fSigma*fSigma));
                    int v0 = ip.getPixel(x,y);
                    ip.putPixelValue(x,y,amplitude*v+v0);
                }//for

        }//end of function drawGaussianDot

        //draw 3D gaussian
        public void draw3DGaussian(ImagePlus imp,int x0,int y0,int z0,int sigma,int amplitude){
            double fSigma = (double)sigma;
            if(fSigma<0.5)fSigma=0.5;//can not be 0!
            int f = 5;
            ImageProcessor ip = imp.getProcessor(); 
            int z1 = z0-sigma*f;
            if(z1<0)z1=0;
            int z2 = z0+sigma*f;
            if(z2>=imp.getStackSize())z2=imp.getStackSize()-1;
            //IJ.log(z0+"="+z1+"::"+z2);

            for(int z = z1;z<=z2;z++){
                imp.setSlice(z);
                //ip.putPixelValue(x0,y0,255);
                for(int y=y0-sigma*f;y<=y0+sigma*f;y++)
                    for(int x=x0-sigma*f;x<=x0+sigma*f;x++){
                        double v =Math.exp((-(x-x0)*(x-x0)-(y-y0)*(y-y0)-(z-z0)*(z-z0))/(2*fSigma*fSigma));
                        int v0 = ip.getPixel(x,y);
                        ip.putPixelValue(x,y,amplitude*v+v0);
                    }//for

            }

        }//end of function drawGaussianDot
        //draw 3D ball
        public void draw3DBall(ImagePlus imp,int x0,int y0,int z0,int sigma,int amplitude){
            double fSigma = (double)sigma;
            if(fSigma<0.5)fSigma=0.5;//can not be 0!
            int f = 5;
            ImageProcessor ip = imp.getProcessor(); 
            ip.setColor(amplitude);
            int z1 = z0-sigma*f;
            if(z1<0)z1=0;
            int z2 = z0+sigma*f;
            if(z2>=imp.getStackSize())z2=imp.getStackSize()-1;
            //IJ.log(z0+"="+z1+"::"+z2);
            double R = sigma*f;

            for(int z = z1;z<=z2;z++){
                imp.setSlice(z);
                double L = z-z0;
                int r = (int)Math.sqrt(R*R-L*L);
                ip.fillOval(x0-r,y0-r,2*r,2*r);
            }

        }//end of function draw3DBall

    }//class ImageTable


    ///////////////////////////////////////////////////////////////
    // This is a mouse wheel listener
    // Mouse wheel can zoom the image
    // Mouse drag can move the display area
    ///////////////////////////////////////////////////////////////
    static class CustomCanvas extends ImageCanvas implements MouseWheelListener { 
        ImageTable itbl;

        // STORM field position (nm) 
        double x0;
        double y0;
        // screen pixel position
        // int screenx0;
        // int screeny0;

        public CustomCanvas(ImagePlus image,ImageTable it) {
            super(image); 
            itbl = it;
            addMouseWheelListener(this); 
        } 
        public void mousePressed(MouseEvent event) {
            int screenx =offScreenX(event.getX());
            int screeny =offScreenY(event.getY());

            x0 = screenx/itbl.zoom+itbl.x_L;
            y0 = screeny/itbl.zoom+itbl.y_U;

            if(event.getButton()==MouseEvent.BUTTON3){//BUTTON3 :Right button
                //screenx0=screenx;
                //screeny0=screeny;
            }


        }

        public void mouseReleased(MouseEvent event) {

            if(event.getButton()==MouseEvent.BUTTON3){//BUTTON3 :Right button
                itbl.createStack();
            }

        }

        public void mouseWheelMoved(MouseWheelEvent event) { 
            synchronized(this) { 
                int wheel = event.getWheelRotation();
                int screenx =offScreenX(event.getX());
                int screeny =offScreenY(event.getY());
                double mouseX = screenx/itbl.zoom+itbl.x_L;
                double mouseY = screeny/itbl.zoom+itbl.y_U;
                double centerX = 0.5*(itbl.x_L+itbl.x_R);
                double centerY = 0.5*(itbl.y_U+itbl.y_D);
                double zoomfactor;
                if (wheel<0) 
                    zoomfactor = 1.1;
                else
                    zoomfactor = 0.9;
                itbl.zoom *=zoomfactor;
                double cx = mouseX*(1.0-1.0/zoomfactor)+centerX/zoomfactor;
                double cy = mouseY*(1.0-1.0/zoomfactor)+centerY/zoomfactor; 
                //IJ.log("cx,cy: ("+event.getX()+","+event.getY()+")");
                //IJ.log("canvas:x="+canvas.getWidth()+" y="+canvas.getHeight());
                itbl.zoomInWindow(itbl.zoom,cx,cy);
            } 
        }

        public void mouseDragged(MouseEvent event) {
            int screenx =offScreenX(event.getX());
            int screeny =offScreenY(event.getY());
            double mouseX = screenx/itbl.zoom+itbl.x_L;
            double mouseY = screeny/itbl.zoom+itbl.y_U;
            double centerX = 0.5*(itbl.x_L+itbl.x_R);
            double centerY = 0.5*(itbl.y_U+itbl.y_D);

            itbl.zoomInWindow(itbl.zoom,x0-mouseX+centerX,y0-mouseY+centerY);
        }

    }//Class CustomCanvas


}//PlugIn STORM_read_table

