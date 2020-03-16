import org.opensourcephysics.numerics.*;


public class SplitOrder {
  double[] Psi,   // Wave Function
             T,   // Kinetic Energy Operator:(dt)       [ Wrapped ]
             V;   // Potential Energy Operator:(dt/2)   [ Natural ]
  FFT2D fft;
  double t, dt, xmin, xmax, ymin, ymax;
  double A;
  int Nx, Ny;

  public SplitOrder(Function f, double xMode, double yMode, double xmin, double max, double ymin, double ymax, int nx, int ny, double xinit, double yinit){
    fft = new FFT2D(nx, ny);     //
    Psi = new double[2*nx*ny];   //  Allocate Memory
    T = new double[2*nx*ny];     //
    V = new double[2*nx*ny];     //
    dt = 1.0/Math.max( fft.getFreqMax(ymin, ymax, ny),  fft.getFreqMax(xmin, xmax, nx));  // minimum allowed timestep
    dt = 0.001;
    t  = 0.0;                    //
    this.xmin = xmin;
    this.xmax = xmax;
    this.ymin = ymin;
    this.ymax = ymax;
    Nx = nx;
    Ny = ny;
    A = Nx*Ny;

    // Psi + V
    double x0 = xmin, dx = (xmax-xmin)/(nx-1);
    double y0 = ymin, dy = (ymax-ymin)/(ny-1);
    // T
    double[] p_x = fft.getWrappedOmegaX(xmin, xmax);
    double[] p_y = fft.getWrappedOmegaY(ymin, ymax);

    for(int iy=0; iy<ny; iy++){
      x0 = xmin;
      int offset = 2*iy*(Nx);
      for(int ix=0; ix<Nx; ix++){
        // Iinitialize psi : (Gaussian Packet)
        double xq = x0 - xinit;
        double yq = y0 + yinit;
        Psi[offset+2*ix] = f.evaluate(xq*xq+yq*yq)*(Math.cos(xMode*xq)*Math.cos(-yMode*yq)
                         -Math.sin(xMode*xq)*Math.sin(-yMode*yq));
        Psi[offset+2*ix+1] = f.evaluate(xq*xq+yq*yq)*(Math.sin(xMode*xq)*Math.cos(-yMode*yq)
                         +Math.cos(xMode*xq)*Math.sin(-yMode*yq));

        // initialize T operator
        double px = p_x[ix], py = p_y[iy];
        T[offset+2*ix] =  Math.cos((px*px+py*py)*dt/2.0);
        T[offset+2*ix+1] = -Math.sin((px*px+py*py)*dt/2.0);

        // initialize V operator
        V[offset+2*ix] =  Math.cos(getV(x0,y0)*dt/2.0);
        V[offset+2*ix+1] = -Math.sin(getV(x0,y0)*dt/2.0);

        x0 += dy;
      }
      y0 += dy;
    }

  }

  double step() {
    double realVT=0, imagVT=0;
    double realP=0, imagP=0;

    for(int i=0; i<Nx; i++){
      int offset = 2*i*Ny;
      for(int j=0; j<Ny; j++){
        realVT = V[offset+2*j];
        imagVT = V[offset+2*j+1];
        realP = Psi[offset+2*j];
        imagP = Psi[offset+2*j+1];

        Psi[offset+2*j] = realVT*realP - imagVT*imagP;
        Psi[offset+2*j+1] = realVT*imagP + imagVT*realP;
      }
    }
    fft.transform(Psi);
    for(int i=0; i<Nx; i++){
      int offset = 2*i*Ny;
      for(int j=0; j<Ny; j++){
        realVT = T[offset+2*j];
        imagVT = T[offset+2*j+1];
        realP = Psi[offset+2*j];
        imagP = Psi[offset+2*j+1];

        Psi[offset+2*j] = realVT*realP - imagVT*imagP;
        Psi[offset+2*j+1] = realVT*imagP + imagVT*realP;
      }
    }
    fft.inverse(Psi);
    for(int i=0; i<Nx; i++){
      int offset = 2*i*Ny;
      for(int j=0; j<Ny; j++){
        realVT = V[offset+2*j];
        imagVT = V[offset+2*j+1];
        realP = Psi[offset+2*j];
        imagP = Psi[offset+2*j+1];

        Psi[offset+2*j] = (realVT*realP - imagVT*imagP)/A;
        Psi[offset+2*j+1] = (realVT*imagP + imagVT*realP)/A;
      }
    }

    return dt;
  }

  public double getV(double x, double y){
    return x*x+y*y;
    /*if(Math.abs(x+1.5)<1.45 && Math.abs(y+1.5)<1.45){
      return 1000.0;
    }
    if(Math.abs(x-1.5)<1.45 && Math.abs(y-1.5)<1.45){
      return 1000.0;
    }
    return 0
    */
  }

}
