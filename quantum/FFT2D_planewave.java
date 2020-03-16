import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.FFT2DFrame;
import org.opensourcephysics.frames.Complex2DFrame;

public class FFT2D_planewave extends AbstractCalculation {
  FFT2DFrame kFrame = new FFT2DFrame("k_x", "k_y", "2D FFT");
  Complex2DFrame xFrame = new Complex2DFrame("x","y","Complex Field");


  public void calculate() {
    int xMode = control.getInt("x mode"),
        yMode = control.getInt("y mode");
    double xmin = control.getDouble("xmin");
    double xmax = control.getDouble("xmax");
    int nx = control.getInt("Nx");
    double ymin = control.getDouble("ymin");
    double ymax = control.getDouble("ymax");
    int ny = control.getInt("Ny");
    // data stored in row-major format
    double[] zdata = new double[2*nx*ny];
    xFrame.setPreferredMinMax(xmin,xmax,ymin,ymax);
    xFrame.resizeGrid(nx,ny);
    xFrame.setAll(zdata);
    double y = 0, yDelta = 2*Math.PI/ny;
    for(int iy=0; iy<ny; iy++){
      // offset to beginning of a row; each row length = nx
      int offset = 2*iy*nx;
      double x = 0, xDelta = 2*Math.PI/nx;
      for(int ix=0; ix<nx; ix++){
        // f(z) = e^(i*xmode*x)e^(i*ymode*y)
        zdata[offset+2*ix] = // real part
           Math.cos(xMode*x)*Math.cos(yMode*y)
          -Math.sin(xMode*x)*Math.sin(yMode*y);
        zdata[offset+2*ix+1] = // imag part
           Math.sin(xMode*x)*Math.cos(yMode*y)
          +Math.cos(xMode*x)*Math.sin(yMode*y);
        x += xDelta;
      }
      y += yDelta;
    }
    xFrame.setAll(zdata);
    xFrame.setVisible(true);
    xFrame.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
    kFrame.doFFT(zdata, nx, xmin, xmax, ymin, ymax);
  }

  public void reset() {
    control.setValue("x mode",0);
    control.setValue("y mode", 1);
    control.setValue("xmin", 0);
    control.setValue("xmax", "2*pi");
    control.setValue("ymin", 0);
    control.setValue("ymax", 1);
    control.setValue("Nx",16);
    control.setValue("Ny",16);
  }

  public static void main(String[] args){
    CalculationControl.createApp(new FFT2D_planewave());
  }
}
