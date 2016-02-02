#ifndef GRADIENT_H
#define GRADIENT_H

#include "itkCovariantVector.h"
#include "itkImage.h"
#include "itkGradientImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"

typedef itk::Image<double,3> ImageType;
typedef itk::CovariantVector< double, 3 > GradientPixelType;
typedef itk::Image< GradientPixelType, 3 > GradientImageType;
typedef itk::GradientImageFilter< ImageType, float, double, GradientImageType > VectorGradientFilterType;
typedef itk::GradientMagnitudeImageFilter< ImageType, ImageType > GradientMFilterType;

class Gradient
{
public:
    Gradient();
    ~Gradient();
    int compute_gradient(ImageType::Pointer image);
    ImageType::Pointer getImageGradient();
    GradientImageType::Pointer getImageVectorGradient();
private:
    //pointeurs itk-> ok
    GradientMFilterType::Pointer m_gradientMagnitudeFilter;
    VectorGradientFilterType::Pointer m_gradientMapFilter;
    
};

#endif // GRADIENT_H
