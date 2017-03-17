package TugasBesarAlgeo;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Scanner;
import static jdk.nashorn.internal.objects.NativeString.toUpperCase;

public class TugasBesarAlgeo
{
    public static void inputmatriks(matriks A)
    {
        int i,j;
        Scanner input = new Scanner(System.in);
        for (i = 0; i < A.NBEff; i++)
        {
            for (j = 0 ; j < A.NKEff ; j++)
            {
                System.out.format("Masukan koefisien ke-%d persamaan ke-%d : ",(j+1),(i+1));
                A.kiri[i][j] = input.nextDouble();
            }
            System.out.format("Masukan konstanta persamaan ke-%d : ",(i+1));
            A.kanan[i] = input.nextDouble();
        }
    }
    public static void cetakmatriks(matriks A)
    {
        int i,j;
        
        for (i = 0; i < A.NBEff; i++)
        {
            for (j = 0 ; j < A.NKEff ; j++)
            {
                System.out.format("%.5f ",A.kiri[i][j]);
                System.out.print("\t");
            }
            System.out.print("| ");
            System.out.format("%.5f\n",A.kanan[i]);
        }
    }
    public static void swap(matriks A,int a,int b)
    {
        int j;
        double temp;
        for (j = 0 ; j < A.NKEff ; j++)
        {
            temp = A.kiri[a][j];
            A.kiri[a][j] = A.kiri[b][j];
            A.kiri[b][j] = temp;
            temp = A.kanan[a];
            A.kanan[a] = A.kanan[b];
            A.kanan[b] = temp;
        }
    }
    public static matriks gauss(matriks A)
    {
        int max,i,j,k,l;
        double bagi;
        for (k = 0, l =0; k < A.NBEff && l< A.NKEff ; k++,l++)
        {
            max = k;
            for (i = (k + 1); i < A.NBEff; i++)
            {
                if ((Math.abs(A.kiri[max][k])) < (Math.abs(A.kiri[i][k])) )
                {
                    max = i;
                }
            }
            swap(A,k,max);
            if (Math.abs(A.kiri[k][k]) <= 0.00000000000000000000000000000000000000000000000000000001)
            {
                matriks temp = A;
                A = temp;
            }
            else
            {
                for (i = (k+1); i< A.NBEff; i++)
                {
                    bagi = A.kiri[i][l] / A.kiri[k][l];
                    for (j = l ;j<A.NKEff;j++)
                    {
                        A.kiri[i][j] = A.kiri[i][j] - A.kiri[k][j] * bagi;
                    }
                    A.kanan[i] = A.kanan[i] - A.kanan[k] * bagi;
                }
            }
        }
        return A;
    }
    public static matriks gaussjordan2(matriks A)
    {
        int max,i,j,k,l;
        double bagi;
        A = gauss(A);
        
        if (jenismatriks(A) == 3)
        {
            for (k = A.NBEff-1, l =A.NKEff-1; k >= 0 && l>= 0 ; k--,l--)
            {

                    for (i = (k-1); i>= 0; i--)
                    {
                        if (A.kiri[k][l] != 0)
                        {
                        bagi = A.kiri[i][l] / A.kiri[k][l];
                        for (j = l ;j<A.NKEff;j++)
                        {
                            A.kiri[i][j] = A.kiri[i][j] - A.kiri[k][j] * bagi;
                        }
                        A.kanan[i] = A.kanan[i] - A.kanan[k] * bagi;
                        }
                    }

            }
            for (k = 0 ; k < A.NBEff ; k++)
            {
                for (l = 0 ; l < A.NKEff ; l++)
                {
                    if ((A.kiri[k][l] != 0) && (k ==l))
                    {
                        A.kanan[k] = A.kanan[k] / A.kiri[k][l];
                        A.kiri[k][l] = A.kiri[k][l]/A.kiri[k][l];
                    }
                }
            }
            
        }
        return A;
    }
    public static matriks gaussjordan(matriks A)
    {
        int i,j;
        A = gauss(A);
        for (i = A.NBEff-1 ; i >=0;i--)
        {
            for (j = A.NKEff-1 ; j>=0;j--)
            {
                if ((j == i ) && (A.kiri[i][j] != 0))
                {
                    A.kanan[i] =  A.kanan[i] / A.kiri[i][j];
                    A.kiri[i][j] = 1;
                } 
                else
                {
                    A.kanan[i] = A.kanan[i] - A.kiri[i][j] * A.kanan[j];
                    A.kiri[i][j] =0;
                }
            }
        }
        return A;
    }
    public static void isiinterpolasi(matriks A, double[] X, double[] Y,int ndata)
    {
        int i,j,k;
        double dummy;
        for (i = 0 ; i< ndata ; i++)
        {
            A.kanan[i] = Y[i];
            for (j = 0 ; j < ndata ;j++)
            {
                if (j == 0 )
                {
                    A.kiri[i][j] = 1; 
                }
                else
                {
                    dummy = 1;
                    for (k = 0 ;k < j ; k++)
                    {
                        dummy = dummy * X[i];
                    }
                    A.kiri[i][j] = dummy;
                }
            }
        }
    }
    public static double hitunginterpolasiX(matriks A,double X)
    {
        double sum,dummy;
        int i,j;
        sum = 0.0;
        A = gaussjordan(A);
        for (i = 0; i<A.NBEff;i++)
        {
            if (i == 0 )
                {
                    dummy = 1; 
                }
                else
                {
                    dummy = 1;
                    for (j = 0 ;j < i ; j++)
                    {
                        dummy = dummy * X;
                    }
                }
            sum = sum + A.kanan[i] * dummy;
        }
        return sum;
    }

    public static matriks bacafile(String ifile)
    {
        matriks A;
        try
        {
            Scanner input = new Scanner(new File(ifile));
            int b = input.nextInt();
            int k = input.nextInt();
            A = new matriks(b,k); 
            while (input.hasNextDouble())
            {
                int i,j;
                for (i = 0; i < A.NBEff; i++)
                {
                    for (j = 0 ; j < A.NKEff ; j++)
                    {
                        A.kiri[i][j] = input.nextDouble();
                    }
                    A.kanan[i] = input.nextDouble();
                }
            }
        }
        catch(Exception e)
        {
           System.out.println("TERJADI KESALAHAN PEMBACAAN FILE EKSTERNAL");
           A = new matriks(0,0);
        }
        return A;
    }
    public static matriks bacainterpolasi(double[] x,double[] y,String ifile,int b)
    {
        try
        {
            Scanner input = new Scanner(new File(ifile));
            b = input.nextInt();
            x = new double[b];
            y = new double[b];
            matriks A = new matriks(b,b);
            while (input.hasNextDouble())
            {
                int i;
                for (i = 0; i < b; i++)
                {
                    x[i] = input.nextDouble();
                }
                for (i = 0; i < b; i++)
                {
                    y[i] = input.nextDouble();
                }
            }
            isiinterpolasi(A,x,y,b);
            return A;
        }
        catch(Exception e)
        {
            System.out.println("TERJADI KESALAHAN PEMBACAAN FILE EKSTERNAL");
            matriks A = new matriks(0,0);
            return A;
        }
    }
    
    public static void cetakmenu()
    {
        boolean inputbenar = false;
        System.out.println("- Aplikasi Aljabar Lanjar pada Metode Numerik -");
        System.out.println("1. Operasi Sistem Persamaan Lanjar");
        System.out.println("2. Interpolasi");
        System.out.print("Masukan Input : ");
        Scanner input = new Scanner(System.in);
        int menu = input.nextInt();
        matriks A = null;
        while (!inputbenar) 
        {
            switch (menu) 
            {
                case 1:
                {
                    while (!(inputbenar)) 
                    {
                        System.out.println("1. Input dari Papan Ketik");
                        System.out.println("2. Input dari File");
                        System.out.print("Masukan Input : ");
                        int menu1 = input.nextInt();
                        switch (menu1)
                        {
                            case 1:
                            {
                                System.out.print("Masukan Jumlah Persamaan : ");
                                int a = input.nextInt();
                                System.out.print("Masukan Jumlah Variabel : ");
                                int b = input.nextInt();
                                A = new matriks(a,b);
                                inputmatriks(A);
                                inputbenar = true;
                                break;
                            }
                            case 2:
                            {
                                System.out.print("Masukan Nama File Input : ");
                                String ifile = input.next();
                                A = bacafile(ifile);
                                inputbenar = true;
                                break;
                            }
                            default:
                            {
                                System.out.println("INPUT SALAH. HARAP INPUT ULANG");
                                inputbenar = false;
                                break;
                            }  
                        }
                    }
                    if ((A.NBEff != 0) && (A.NKEff != 0))
                    {
                        System.out.println("Matriks dari Persamaan yang Diinput : ");
                        cetakmatriks(A);
                        inputbenar = false;
                        while (!(inputbenar))
                        {
                            System.out.println("Operasi SPL yang Diinginkan : ");
                            System.out.println("1. Eliminasi Gauss");
                            System.out.println("2. Eliminasi Gauss-Jordan");
                            System.out.print("Masukan Input : ");
                            int menu2 = input.nextInt();
                            switch (menu2) 
                            {
                                case 1:
                                {
                                    matriks B = A;
                                    B = gauss(B);
                                    cetakmatriks(B);
                                    inputbenar = true;
                                    break;
                                }
                                case 2:
                                {
                                    matriks B = A;
                                    B = gaussjordan2(B);
                                    cetakmatriks(B);
                                    inputbenar = true;
                                    break;
                                }
                                default:
                                {
                                    System.out.println("INPUT SALAH. HARAP INPUT ULANG");
                                    inputbenar = false;
                                    break;
                                }
                            }
                        }
                        cetakhasil(A);
                        System.out.println("Apakah Mau Disimpan Ke File Eksternal ?");
                        System.out.print("Masukan Input : ");
                        String menu3 = input.next();
                        if ("YA".equals(toUpperCase(menu3)))
                        {
                            tulismatriksfile(A);
                        }
                    }
                    break;
                }
                
                case 2:
                {
                    int a = 0;
                    double[] x = null;
                    double [] y = null;
                    while (!(inputbenar)) 
                    {
                        System.out.println("1. Input dari Papan Ketik");
                        System.out.println("2. Input dari File");
                        System.out.print("Masukan Input : ");
                        int menu1 = input.nextInt();
                        switch (menu1)
                        {
                            case 1:
                            {
                                System.out.print("Masukan Jumlah Data : ");
                                a = input.nextInt();
                                x = new double[a];
                                y = new double[a];
                                int i;
                                for (i = 0;i< a ; i++)
                                {
                                    System.out.print("Masukan Nilai X : ");
                                    x[i] = input.nextDouble();
                                    System.out.print("Masukan Nilai F(X) : ");
                                    y[i] = input.nextDouble();
                                }
                                A = new matriks(a,a);
                                isiinterpolasi(A,x,y,a);
                                inputbenar = true;
                                break;
                            }
                            case 2:
                            {
                                System.out.print("Masukan Nama File Input : ");
                                String ifile = input.next();
                                A = bacainterpolasi(x,y,ifile,a);
                                inputbenar = true;
                                break;
                            }
                            default:
                            {
                                System.out.println("INPUT SALAH. HARAP INPUT ULANG");
                                inputbenar = false;
                                break;
                            }  
                        }
                    }
                    if ((A.NBEff != 0) && (A.NKEff != 0))
                    {
                        System.out.println("Matriks dari Data Interpolasi : ");
                        cetakmatriks(A);
                        System.out.println();
                        cetakmatriks(gaussjordan(A));
                        System.out.println("Persamaan polinomnya : ");
                        System.out.println(hasilpolinom(A));
                        System.out.print("Masukan Nilai X yang Ingin Hitung dengan Data Interpolasi : ");
                        double nilaiX = input.nextDouble();
                        double hasil = hitunginterpolasiX(A,nilaiX);
                        System.out.format("F(%.5f) = %.5f\n",nilaiX,hasil);
                        System.out.println("Apakah Mau Disimpan Ke File Eksternal ?");
                        System.out.print("Masukan Input : ");
                        String menu3 = input.next();
                        if ("YA".equals(toUpperCase(menu3)))
                        {
                            tulishasilpolinom(hasilpolinom(A));
                        }
                    }
                    break;
                }
                default:
                {
                    System.out.println("INPUT SALAH. HARAP INPUT ULANG");
                    inputbenar = false;
                    break ;
                }
            }
        }
    }
    public static String hasilpolinom(matriks A)
    {
        int i, j;
        String z = "";
        z = z + "F(X) = ";
        for (i=0; i<A.NBEff; i++) 
        {
            if (A.kanan[i] != 0)
            {
                if (A.kanan[i] != 1)
                {
                    if (A.kanan[i] > 0) 
                    {
                        if (i != 0)
                        {
                            z = z + " + ";
                        }
                        z = z + String.valueOf(A.kanan[i]);
                    }
                    else if (A.kanan[i] < 0)
                    {
                        z = z + " - ";
                        z = z + String.valueOf(-1 * A.kanan[i]);
                    } 
                }
                else if ((A.kanan[i]==1) && (i ==0))
                {
                    z = z + "1";
                }
                else if ((A.kanan[i]==1) && (i !=0))
                {
                    z = z + " + ";
                }
                switch (i) {
                    case 1:
                        z = z + "X";
                        break;
                    case 0:
                        break;
                    default:
                        z = z + "X^";
                        z = z + String.valueOf(i);
                        break;
                }
            }
        }
        System.out.println();
        return z;
    }
    public static void tulismatriksfile(matriks A)
    {
        try {
            
            FileWriter fileWriter = new FileWriter("outputmatriks.txt");
            
            
            BufferedWriter bw = new BufferedWriter(fileWriter);
            int i,j;
            File output = new File ("outputmatriks.txt");
            if (!output.exists())
                output.createNewFile();
            bw.write("Matriks :");
            bw.newLine();
            for (i = 0;i< A.NBEff;i++)
            {
                for (j =0;j< A.NKEff; j++)
                {
                    bw.write(Double.toString(A.kiri[i][j])) ;
                    bw.write("\t");
                }
                bw.write("\t|\t");
                bw.write(Double.toString(A.kanan[i]));
                bw.newLine();
            }
            bw.write("Hasil setelah di Gauss Jordan :");
            bw.newLine();
            matriks B = gaussjordan2(A);
            for (i = 0;i< A.NBEff;i++)
            {
                for (j =0;j< A.NKEff; j++)
                {
                    bw.write(Double.toString(B.kiri[i][j])) ;
                    bw.write("\t");
                }
                bw.write("\t|\t");
                bw.write(Double.toString(B.kanan[i]));
                bw.newLine();
            }
            bw.close();
        }
        catch(Exception ex) {
            System.out.println("TERJADI KESALAHAN PADA PENULISAN FILE");
        }
    }
    public static void tulishasilpolinom(String z)
    {
        try {
            
            FileWriter fileWriter = new FileWriter("outputinterpolasi.txt");
            
            
            BufferedWriter bw = new BufferedWriter(fileWriter);
            int i,j;
            File output = new File ("outputinterpolasi.txt");
            if (!output.exists())
                output.createNewFile();
            bw.write("Polinom Interpolasi");
            bw.newLine();
            bw.write(z);
            bw.close();
        }
        catch(Exception ex) {
            System.out.println("TERJADI KESALAHAN PADA PENULISAN FILE");
        }
    }
    public static void cetakhasil(matriks A)
    {
        A=gauss(A);
        switch(jenismatriks(A))
        {
            case 1 :
            {
                System.out.println("Sistem Persamaan Lanjar Memiliki Banyak Solusi");
                break;
            }
            case 2 :                
            {
                System.out.println("Sistem Persamaan Lanjar Tidak Memiliki Solusi");
                break;
            }
            default :
            {
                A = gaussjordan2(A);
                System.out.println("Sistem Persamaan Lanjar Memiliki Solusi Unik :");
                int i;
                for (i = 0; i< A.NBEff;i++)
                {
                    System.out.format("X[%d] = %.5f \n",(i+1),A.kanan[i]);
                }
                System.out.println();
                break;
            }
        }
    }
    public static int jenismatriks(matriks A)
    {
        int j,j0;
        j0 = 0;
        j=0;
        while((A.kiri[A.NBEff-1][j] == 0) && (j< A.NKEff))
        {
            if (A.kiri[A.NBEff-1][j] == 0) 
            {
                j0++;
            }
            j++;
        }
        if (j0 == (A.NKEff))
        {
            if (A.kanan[A.NBEff-1] == 0)
            {
                return 1;
            }
            else
            {
                return 2;
            }
        }
        else if (j0 == (A.NKEff-1))
        {
            j0=0;
            j=0;
            while((A.kiri[A.NBEff-1][j] == 0) && (j< A.NKEff))
            {
                if (A.kiri[A.NBEff-2][j] == 0) 
                {
                    j0++;
                }
                j++;
            }
            if (j0 == (A.NKEff-1))
            {
                return 2;
            }
            else if (j0 == (A.NKEff-2))
            {
                return 3;
            }
            else
            {
                return 1;
                        
            }
        }
        else
        {
            return 1;
        }
    }
    public static void main(String[] args) 
    {
        cetakmenu();
    }
    
}
