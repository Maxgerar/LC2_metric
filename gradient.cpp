#include "gradient.h"

Gradient::Gradient()
{
    m_gradientMagnitudeFilter = GradientMFilterType::New() ;
    m_gradientMapFilter = VectorGradientFilterType::New() ;
    
}

Gradient::~Gradient()
{
    
}

int Gradient::compute_gradient(ImageType::Pointer image)
{
    m_gradientMagnitudeFilter->SetInput(image);
    try {
        m_gradientMagnitudeFilter->Update();
    } catch( itk::ExceptionObject & e ) {
        std::cerr << "Exception caught while updating gradientMagnitudeFilter " << std::endl;
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    
    // vector gradient -> contains vector gradient for each pixel
    m_gradientMapFilter->SetInput( image );
    try {
        m_gradientMapFilter->Update();
    } catch( itk::ExceptionObject & e ) {
        std::cerr << "Exception caught while updating gradientMapFilter " << std::endl;
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    return 0;
}

ImageType::Pointer Gradient::getImageGradient()
{
    return m_gradientMagnitudeFilter->GetOutput();
}

GradientImageType::Pointer Gradient::getImageVectorGradient()
{
    return m_gradientMapFilter->GetOutput();
}

