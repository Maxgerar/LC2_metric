#ifndef LC2_H
#define LC2_H

#include "itkimage.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageToImageMetric.h"

#include <vector>

typedef itk::Image<double,3> ImageType;
typedef itk::ImageToImageMetric<ImageType, ImageType> MetricType;
//objet matrice et vecteur pour le calcul des parametre a,b et c
typedef itk::Matrix<double, 343,3> MType;
typedef itk::Matrix<double,3,343> MTType;
typedef itk::Matrix<double,3,3> M3Type;
typedef itk::Vector<double, 343> VType;
typedef itk::Vector<double,3> PType;
//s = 3 -> m = 343 voxels dans un voisinage !
//voisinage
typedef itk::NeighborhoodIterator<ImageType> NIteratorType;
typedef itk::ImageRegionConstIterator<ImageType> ImageConstIteratorType;

using namespace std;

//mes fichiers
#include "gradient.h"

class Lc2: public MetricType { //:public MetricType
public:
  
    typedef Lc2 Self;
    typedef ImageToImageMetric< ImageType, ImageType > Superclass;
    typedef itk::SmartPointer<Self>                             Pointer;
    typedef itk::SmartPointer<const Self>                    ConstPointer;
    
    //itkNewMacro(Self);
    itkTypeMacro(Lc2, Object);
    
    typedef typename Superclass::TransformParametersType TransformParametersType;
    typedef typename Superclass::MeasureType             MeasureType;

    
    
    MeasureType GetValue(const TransformParametersType);
    double getlc2();
    void compute_param();
    
    void setImages(ImageType::Pointer image_IRM,ImageType::Pointer image_US);
    
protected:
    //constructeur
    Lc2();
	
    
private:
	double lc2final;
    double variancesum;
    double lc2varsum;
    ImageType::Pointer im_US;
    ImageType::Pointer im_IRM;
    ImageType::Pointer im_grad_IRM;
    vector<vector<double>> parameters;
    vector<double> variances;
    vector<double> LC2_i;
    
};
#endif