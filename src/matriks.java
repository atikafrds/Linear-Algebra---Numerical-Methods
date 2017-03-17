package TugasBesarAlgeo;
public class matriks 
{
    public int NBEff;
    public int NKEff;
    public double[][] kiri;
    public double[] kanan;
    public matriks(int a, int b) 
    {
        NBEff = a;
        NKEff = b;
        this.kiri = new double[NBEff][NKEff];
        this.kanan = new double[NBEff];
    }
}
