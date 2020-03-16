import org.opensourcephysics.controls.*;
import org.opensourcephysics.display.*;
import org.opensourcephysics.frames.*;
import org.opensourcephysics.numerics.*;

public class FFT2DCalculationApp extends AbstractCalculation {
  //FFT2DFrame kFrame = new FFT2DFrame("k_x", "k_y", "2D FFT");
  Complex2DFrame xFrame = new Complex2DFrame("x","y","Complex Field");
  Complex2DFrame pFrame = new Complex2DFrame("p_x","p_y","Momentum Space");

  public void calculate() {
    xFrame.clearData();
    pFrame.clearData();
    String fStr = control.getString("f(|r|)");
    Function f = null;
    try{
      f = new ParsedFunction(fStr, "r");
    } catch(ParserException ex){
      control.println("Error parsing function string: "+fStr);
      return;
    }
    int xMode = control.getInt("x mode"),
        yMode = control.getInt("y mode");
    double xmin = control.getDouble("xmin");
    double xmax = control.getDouble("xmax");
    int nx = control.getInt("Nx");
    double ymin = control.getDouble("ymin");
    double ymax = control.getDouble("ymax");
    int ny = control.getInt("Ny");

   FFT2D fft = new FFT2D(nx, ny);
   double px_min = fft.getFreqMin(xmin, xmax, nx)*Math.PI*2.0;
   double px_max = fft.getFreqMax(xmin, xmax, nx)*Math.PI*2.0;
   double py_min = fft.getFreqMin(ymin, ymax, ny)*Math.PI*2.0;
   double py_max = fft.getFreqMax(ymin, ymax, ny)*Math.PI*2.0;

    // data stored in row-major format
    double[] zdata = new double[2*nx*ny];
    xFrame.setPreferredMinMax(xmin,xmax,ymin,ymax);
    xFrame.resizeGrid(nx,ny);
    xFrame.setAll(zdata);
    pFrame.setPreferredMinMax(px_min,px_max,py_min,py_max);
    pFrame.resizeGrid(nx,ny);
    pFrame.setAll(zdata);
    for(int iy=0; iy<ny; iy++){
      // offset to beginning of a row; each row length = nx
      int offset = 2*iy*nx;
      double y = xFrame.indexToY(iy);
      for(int ix=0; ix<nx; ix++){
        double x = xFrame.indexToX(ix);
        // f(z) = e^(i*xmode*x)e^(i*ymode*y)
        zdata[offset+2*ix] = // real part
           f.evaluate(x*x+y*y)*(Math.cos(xMode*x)*Math.cos(yMode*y)
          -Math.sin(xMode*x)*Math.sin(yMode*y));
        zdata[offset+2*ix+1] = // imag part
           f.evaluate(x*x+y*y)*(Math.sin(xMode*x)*Math.cos(yMode*y)
          +Math.cos(xMode*x)*Math.sin(yMode*y));
      }

    }

    //new double[2*nx*ny];
    /*
    int offX = (int) (nx*xmin/(xmax-xmin));
    offX = Math.abs(offX);
    int offY = (int) (ny*ymin/(ymax-ymin));
    offY = Math.abs(offY);
    for(int j = 0; j<ny; j++) {
      int jj = (offY+j)%ny;

      int offset = 2*j*nx;
      int offset2 = 2*jj*nx;

      for(int i = 0; i<nx; i++) {
        int ii = (offX+i)%nx;

        Zdata[offset+2*ii] = zdata[offset2+2*i];
        Zdata[offset+2*ii+1] = zdata[offset2+2*i+1];
      }

    }
    */

    //kFrame.setDomainType(4);
    //kFrame.doFFT(zdata, nx, xmin, xmax, ymin, ymax);

    double[] p_x = fft.getNaturalOmegaX(xmin, xmax);
    double[] p_y = fft.getNaturalOmegaY(ymin, ymax);
    double[] p_x2 = fft.getWrappedOmegaX(xmin, xmax);
    double[] p_y2 = fft.getWrappedOmegaY(ymin, ymax);

    double dt = 1.0/Math.max( fft.getFreqMax(ymin, ymax, ny),  fft.getFreqMax(xmin, xmax, nx));// 0.01;
    // T[i][j] = e^(i*dt*p^2/2);
    double[] T = new double[2*nx*ny];
    double[] T2 = new double[2*nx*ny];
    for(int iy=0; iy<ny; iy++){
      int offset = 2*iy*nx;
      double py = p_y[iy];
      double py2 = p_y2[iy];
      for(int ix=0; ix<nx; ix++){
        double px = p_x[ix];
        double px2 = p_x2[ix];
        // real part
        T[offset+2*ix]   = Math.cos((px*px+py*py)*dt/2.0);
        // imag part
        T[offset+2*ix+1] = -Math.sin((px*px+py*py)*dt/2.0);
        // real part
        T2[offset+2*ix]   = Math.cos((px2*px2+py2*py2)*dt/2.0);
        // imag part
        T2[offset+2*ix+1] = -Math.sin((px2*px2+py2*py2)*dt/2.0);

      }
    }

    double[] Zdata = zdata.clone();


    fft.transform(Zdata);


    fft.inverse(Zdata);

    pFrame.setAll(Zdata);

    fft.transform(zdata);
    fft.backtransform(zdata);


    xFrame.setAll(zdata);

    xFrame.setVisible(true);
    xFrame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
    pFrame.setVisible(true);
    pFrame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);

  }

  public void reset() {
    control.setValue("f(|r|)", "exp(-r*r)");
    control.setValue("x mode",0);
    control.setValue("y mode", 5);
    control.setValue("xmin", "-pi");
    control.setValue("xmax", "pi");
    control.setValue("ymin", "-pi");
    control.setValue("ymax", "pi");
    control.setValue("Nx",32);
    control.setValue("Ny",32);
  }

  public static void main(String[] args){
    CalculationControl.createApp(new FFT2DCalculationApp());
  }
}
