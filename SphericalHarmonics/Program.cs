using System;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.Serialization.Formatters.Binary;
using MathNet.Numerics;

namespace SphericalHarmonics
{
    class Program
    {
        // store precalculated associated legendre polynomial values: P_lm(cos(theta))
        static double[][][] legendres;

        static void Main(string[] args)
        {
            // width and height of resulting images (height must be same as pre-generated Octave data
            int width = 6000;
            int height = 3000;
            // if legendres are calculated between 0 and 49 this should be 50, how far up to generate images
            int max_l = 50;
            // where the pre-generated data is stored
            string csvPath = "C:\\Users\\bouellet\\octave scripts";
            legendres = readLegendres(max_l, csvPath);
            // generate rectilinear projection of example spherical harmonic temperature oscillations (brute force)
            for (int l = 1; l <= max_l; l++)
            {
                for (int m = -l; m <= l; m++)
                {
                    double[,] data = new double[width, height];
                    for (int x = 0; x < width; x++)
                    {
                        for (int y = 0; y < height; y++)
                        {
                            // theta (vertical and varies between 0 -> pi)
                            // phi (horizontal and varies between 0 -> 2*pi)
                            // x,y for C# bitmap is upper left corner
                            double theta = Math.PI * y / height;
                            double phi = 2 * Math.PI * x / width;
                            // since theta is only used as the argument to the associated legendre polynomial
                            // and the Octave script has already encoded the cos(theta) into 0->height samples
                            // we are using y directly to index the Octave data
                            data[x,y] = ValueFromPosition(phi, y, l, m);
                        }
                    }
                    // get min and max of all data to scale colors appropriately
                    double min = data.Cast<double>().Min();
                    double max = data.Cast<double>().Max();
                    // set each pixel in a new image to be equal to the calculated color, mapped between blue and red (low to high)
                    Bitmap finalBmp = new Bitmap(width, height);
                    for (int x = 0; x < width; x++)
                    {
                        for (int y = 0; y < height; y++)
                        {
                            finalBmp.SetPixel(x, y, MapColor(data[x,y], min, max));
                        }
                    }
                    // output the final image
                    finalBmp.Save("sphericalharmonic-" + l + "_" + m + ".png", System.Drawing.Imaging.ImageFormat.Png);
                }
            }
        }

        static double[][][] readLegendres(int max_l, string csvPath)
        {
            // since C# and Math.NET don't have associated legendre polynomials I went ahead and
            // precomputed unnormalized associated legendre polynomials in Octave for l between 0 and 49
            // using cos with a sampling equal to the height of the desired image, script below:
            /*
             *  for i = 0:49 
             *      csvwrite(strcat(strcat("legendres-",num2str(i)),".csv"), legendre (i, cos(0:pi/3000:pi*2999/3000)));
             *  endfor
             */
            // serialize out to a binary file to avoid reparsing csv every time this runs
            if (File.Exists("legendres.bin"))
            {
                using (Stream stream = new FileStream("legendres.bin", FileMode.Open))
                {
                    var formatter = new BinaryFormatter();
                    return (double[][][])formatter.Deserialize(stream);
                }
            }
            double[][][] legendres = new double[max_l][][];
            // file is generated for each l value and saved in filename
            foreach (var file in Directory.EnumerateFiles(csvPath, "legendres-*.csv"))
            {
                int l = int.Parse(Path.GetFileNameWithoutExtension(file).Replace("legendres-", ""));
                legendres[l] = new double[l+1][];
                int m = 0;
                // parse in CSV data as doubles
                var parser = new Microsoft.VisualBasic.FileIO.TextFieldParser(file);
                parser.TextFieldType = Microsoft.VisualBasic.FileIO.FieldType.Delimited;
                parser.SetDelimiters(new string[] { "," });
                while (!parser.EndOfData)
                {
                    // parse entire row of values into a double array for this l,m combination
                    // each subsequent row in the octave results represents the next m
                    string[] row = parser.ReadFields();
                    legendres[l][m] = row.Select(i => double.Parse(i, NumberStyles.Float)).ToArray();                    
                    m++;
                }
            }
            using (Stream stream = new FileStream("legendres.bin", FileMode.Create))
            {
                var formatter = new BinaryFormatter();
                formatter.Serialize(stream, legendres);
            }
            return legendres;
        }
        
        static double ValueFromPosition(double phi, int y, int l, int m)
        {
            // http://find.spa.umn.edu/~pryke/logbook/20000922/
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
            // http://limbicsoft.com/volker/prosem_paper.pdf
            // https://www.physicsforums.com/threads/can-someone-explain-angular-power-spectrum.309483/
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Negative_m_and.2For_negative_.E2.84.93
            // l >= 0
            // m >= 0
            // m <= l
            // for a single component of harmonic decomposition:
            // T = a_l,m * Y = a_l,m * n_l,m * P_l,m(cos(theta)) * e^(i * m * phi)
            // Re(T) = Re(a_l,m * Y) = n_l,m * P_l,m(cos(theta)) * (Real(a_l,m) * cos(m * phi) - Imag(a_l,m) * sin(m * phi))
            // if m is negative P_l,-m = (-1)^m * (l - m)! / (l + m)! * P_l,m
            // a_l,m is removed from the equation for display
            var result = CalcN(l, m) * legendres[l][Math.Abs(m)][y] * (Math.Cos(m * phi) - Math.Sin(m * phi));
            if (m < 0)
            {
                result *= Math.Pow(-1, m) * SpecialFunctions.Factorial(l - m) / SpecialFunctions.Factorial(l + m);
            }
            return result;
        }

        static double CalcN(int l, int m)
        {
            // n_l,m = sqrt((2 * l + 1) / (4 * pi) * (l - m)! / (l + m)!)
            return Math.Sqrt((2 * l + 1) / (4 * Math.PI) * SpecialFunctions.Factorial(l - m) / SpecialFunctions.Factorial(l + m));
        }

        // map a value to a color
        static Color MapColor(double value, double blue_value, double red_value)
        {
            // Convert into a value between 0 and 1023.
            int int_value = (int)(1023 * (value - blue_value) / (red_value - blue_value));
            if (red_value == blue_value)
            {
                int_value = 1024 / 2;
            }
            // Map different color bands.
            if (int_value < 256)
            {
                // blue to aqua: (0, 0, 255) to (0, 255, 255)
                return Color.FromArgb(0, int_value, 255);
            }
            else if (int_value < 512)
            {
                // aqua to green: (0, 255, 255) to (0, 255, 0)
                int_value -= 256;
                return Color.FromArgb(0, 255, 255 - int_value);
            }
            else if (int_value < 768)
            {
                // green to yellow: (0, 255, 0) to (255, 255, 0)
                int_value -= 512;
                return Color.FromArgb(int_value, 255, 0);
            }
            else
            {
                // yellow to red: (255, 255, 0) to (255, 0, 0)
                int_value -= 768;
                return Color.FromArgb(255, 255 - int_value, 0);
            }
        }
    }
}
