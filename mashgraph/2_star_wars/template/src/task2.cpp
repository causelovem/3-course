#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cmath>

#include "classifier.h"
#include "EasyBMP.h"
#include "linear.h"
#include "argvparser.h"
#include "matrix.h"

#define PI 3.1415926535897932384626433832795
#define EPS 0.001
#define SEG_NUM 16
#define BLOCK_SIZE 5

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;

using CommandLineProcessing::ArgvParser;

typedef vector<pair<BMP*, int> > TDataSet;
typedef vector<pair<string, int> > TFileList;
typedef vector<pair<vector<float>, int> > TFeatures;

class UnaryMatrixOp /*+++++*/
{
    public:
        Matrix <float> kernel;
        size_t vert_radius;
        size_t hor_radius;

        UnaryMatrixOp(Matrix<float> kern):
        kernel{kern}, vert_radius{(kern.n_rows - 1) / 2}, hor_radius{(kern.n_cols - 1) / 2}
        {};

        float operator () (const Matrix <float> &img)
        {
            float s = 0;

            for (size_t i = 0; i < img.n_rows; i++)
                for (size_t j = 0; j < img.n_cols; j++) 
                    s += img(i,j) * kernel(i,j);

            /*if (s < 0)
                s = 0;
            if (s > 255)
                s = 255;*/

            return s;
        }
};

// Load list of files and its labels from 'data_file' and
// stores it in 'file_list'
void LoadFileList(const string& data_file, TFileList* file_list)
{
    ifstream stream(data_file.c_str());

    string filename;
    int label;
    
    int char_idx = data_file.size() - 1;
    for (; char_idx >= 0; --char_idx)
        if (data_file[char_idx] == '/' || data_file[char_idx] == '\\')
            break;
    string data_path = data_file.substr(0,char_idx+1);
    
    while(!stream.eof() && !stream.fail()) {
        stream >> filename >> label;
        if (filename.size())
            file_list->push_back(make_pair(data_path + filename, label));
    }

    stream.close();
}

// Load images by list of files 'file_list' and store them in 'data_set'
void LoadImages(const TFileList& file_list, TDataSet* data_set)
{
    for (size_t img_idx = 0; img_idx < file_list.size(); ++img_idx) {
            // Create image
        BMP* image = new BMP();
            // Read image from file
        image->ReadFromFile(file_list[img_idx].first.c_str());
            // Add image and it's label to dataset
        data_set->push_back(make_pair(image, file_list[img_idx].second));
    }
}

// Save result of prediction to file
void SavePredictions(const TFileList& file_list,
                     const TLabels& labels, 
                     const string& prediction_file)
{
        // Check that list of files and list of labels has equal size 
    assert(file_list.size() == labels.size());
        // Open 'prediction_file' for writing
    ofstream stream(prediction_file.c_str());

        // Write file names and labels to stream
    for (size_t image_idx = 0; image_idx < file_list.size(); ++image_idx)
        stream << file_list[image_idx].first << " " << labels[image_idx] << endl;
    stream.close();
}

// Exatract features from dataset.
// You should implement this function by yourself =)
void ExtractFeatures(const TDataSet& data_set, TFeatures* features)
{
    std::vector<pair<float, float>> segments;
    segments.resize(SEG_NUM);

    int denom = -SEG_NUM / 2, count = 0;

    for (int i = denom; i < 0; i++)
        std::get<0>(segments[count++]) = i * PI / denom;

    denom *= -1;
    for (int i = 1; i < denom + 1; i++)
        std::get<0>(segments[count++]) = i * PI / denom;

    Matrix <float> sobel_x(1,3);
    Matrix <float> sobel_y(3,1);

    sobel_x(0, 0) = -1.0;
    sobel_x(0, 1) = 0;
    sobel_x(0, 2) = 1.0;

    sobel_x(0, 0) = -1.0;
    sobel_x(1, 0) = 0;
    sobel_x(2, 0) = 1.0;

    for (size_t image_idx = 0; image_idx < data_set.size(); image_idx++)
    {
        BMP *img;
        img = std::get<0>(data_set[image_idx]);
        Matrix <float> bright(img->TellHeight(), img->TellWidth());

        for (int i = 0; i < img->TellHeight(); i++)
            for (int j = 0; j < img->TellWidth(); j++)
                bright(i,j) = 0.299 * img->GetPixel(i,j).Red + 0.587 * img->GetPixel(i,j).Green + 0.114 * img->GetPixel(i,j).Blue;

        Matrix <float> Ix(img->TellHeight(), img->TellWidth());
        Matrix <float> Iy(img->TellHeight(), img->TellWidth());

        UnaryMatrixOp opx(sobel_x);
        UnaryMatrixOp opy(sobel_y);
        Ix = bright.unary_map(opx);
        Iy = bright.unary_map(opy);

        Matrix <float> modG(Ix.n_rows, Ix.n_cols);
        Matrix <float> theta(Ix.n_rows, Ix.n_cols);

        for (size_t i = 0; i < Ix.n_rows; i++)
            for (size_t j = 0; j < Ix.n_cols; j++)
            {
                modG(i, j) = sqrt (Ix(i, j) * Ix(i, j) + Iy(i, j) * Iy(i, j));
                theta(i, j) = (-PI / 2) + (PI / 4) * round((4 / PI) * atan2(Iy(i, j), Ix(i, j)));
            }

        


        
        /*// PLACE YOUR CODE HERE
        // Remove this sample code and place your feature extraction code here
        vector<float> one_image_features;
        one_image_features.push_back(1.0);
        features->push_back(make_pair(one_image_features, 1));
        // End of sample code*/
    }
}

// Clear dataset structure
void ClearDataset(TDataSet* data_set)
{
        // Delete all images from dataset
    for (size_t image_idx = 0; image_idx < data_set->size(); ++image_idx)
        delete (*data_set)[image_idx].first;
        // Clear dataset
    data_set->clear();
}

// Train SVM classifier using data from 'data_file' and save trained model
// to 'model_file'
void TrainClassifier(const string& data_file, const string& model_file)
{
        // List of image file names and its labels
    TFileList file_list;
        // Structure of images and its labels
    TDataSet data_set;
        // Structure of features of images and its labels
    TFeatures features;
        // Model which would be trained
    TModel model;
        // Parameters of classifier
    TClassifierParams params;
    
        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);

        // PLACE YOUR CODE HERE
        // You can change parameters of classifier here
    params.C = 0.01;
    TClassifier classifier(params);
        // Train classifier
    classifier.Train(features, &model);
        // Save model to file
    model.Save(model_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

// Predict data from 'data_file' using model from 'model_file' and
// save predictions to 'prediction_file'
void PredictData(const string& data_file,
                 const string& model_file,
                 const string& prediction_file)
{
        // List of image file names and its labels
    TFileList file_list;
        // Structure of images and its labels
    TDataSet data_set;
        // Structure of features of images and its labels
    TFeatures features;
        // List of image labels
    TLabels labels;

        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);

        // Classifier 
    TClassifier classifier = TClassifier(TClassifierParams());
        // Trained model
    TModel model;
        // Load model from file
    model.Load(model_file);
        // Predict images by its features using 'model' and store predictions
        // to 'labels'
    classifier.Predict(features, model, &labels);

        // Save predictions
    SavePredictions(file_list, labels, prediction_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

int main(int argc, char** argv)
{
    // Command line options parser
    ArgvParser cmd;
        // Description of program
    cmd.setIntroductoryDescription("Machine graphics course, task 2. CMC MSU, 2014.");
        // Add help option
    cmd.setHelpOption("h", "help", "Print this help message");
        // Add other options
    cmd.defineOption("data_set", "File with dataset",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("model", "Path to file to save or load model",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("predicted_labels", "Path to file to save prediction results",
        ArgvParser::OptionRequiresValue);
    cmd.defineOption("train", "Train classifier");
    cmd.defineOption("predict", "Predict dataset");
        
        // Add options aliases
    cmd.defineOptionAlternative("data_set", "d");
    cmd.defineOptionAlternative("model", "m");
    cmd.defineOptionAlternative("predicted_labels", "l");
    cmd.defineOptionAlternative("train", "t");
    cmd.defineOptionAlternative("predict", "p");

        // Parse options
    int result = cmd.parse(argc, argv);

        // Check for errors or help option
    if (result) {
        cout << cmd.parseErrorDescription(result) << endl;
        return result;
    }

        // Get values 
    string data_file = cmd.optionValue("data_set");
    string model_file = cmd.optionValue("model");
    bool train = cmd.foundOption("train");
    bool predict = cmd.foundOption("predict");

        // If we need to train classifier
    if (train)
        TrainClassifier(data_file, model_file);
        // If we need to predict data
    if (predict) {
            // You must declare file to save images
        if (!cmd.foundOption("predicted_labels")) {
            cerr << "Error! Option --predicted_labels not found!" << endl;
            return 1;
        }
            // File to save predictions
        string prediction_file = cmd.optionValue("predicted_labels");
            // Predict data
        PredictData(data_file, model_file, prediction_file);
    }
}