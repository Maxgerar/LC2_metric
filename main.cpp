//
//  main.cpp
//  LC2_M
// Calcul de la metrique LC2
//  Created by Maxime Gérard on 17/12/15.
//  Copyright © 2015 Maxime Gérard. All rights reserved.
//

#include <iostream>
#include <string>
#include "gradient.h"
#include "Lc2.h"


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegistrationMethod.h"
#include "itkAmoebaOptimizer.h"
#include "itkNiftiImageIO.h"
#include "itkPNGImageIO.h"
#include "itkMetaImageIO.h"
#include "itkImageIterator.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRigid3Dtransform.h"
#include "itkAmoebaOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkEuler3DTransform.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include"itkLC2ImageToImageMetric.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkShrinkImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkCommand.h"




//images
typedef itk::Image<double,3> ImageType;
typedef itk::Image<double,2> Image2DType;

//IO
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileReader<Image2DType> Reader2DType;
typedef itk::ImageFileWriter<ImageType> WriterType;

//iteration over images
typedef itk::ImageRegionIterator<ImageType> IteratorType;

//recalage
typedef itk::AmoebaOptimizer OptimizerType;
//typedef itk::LevenbergMarquardtOptimizer OptimizerType;
typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;
typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
typedef itk::Euler3DTransform<double> TransformType;
typedef itk::LC2ImageToImageMetric<ImageType, ImageType> LC2MetricType;
typedef itk::TranslationTransform<double,3> TranslationType;


//pour le mapping de l'image mobile registree ac tsf determinee par registration framework
typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;

//pour la mise a l'echelle de la bande passante des intensites
typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> RescaleFilterType;
typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxCalculatorType;
typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkFilterType;


using namespace std;

//classe pour suivre le recalage
class CommandIterationUpdate : public itk::Command
{
public :
    typedef CommandIterationUpdate Self;
    typedef itk::Command SuperClass;
    typedef itk::SmartPointer<Self> Pointer;
    itkNewMacro(Self);
    
protected:
    CommandIterationUpdate()
    {
        m_IterationNumber =0;
    }
public:
    typedef itk::AmoebaOptimizer OptimizerType;
    typedef const OptimizerType * OptimizerPointer;
    
    void Execute(itk::Object *caller, const itk::EventObject &event) ITK_OVERRIDE
    {
        Execute( (const itk::Object *)caller, event);
    }
    
    void Execute(const itk::Object * object,
                 const itk::EventObject & event) ITK_OVERRIDE
    {
        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
        if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
            return;
        }
        std::cout << m_IterationNumber++ << "   "<<endl;
        std::cout << optimizer->GetCachedValue() << "   "<<endl;
        std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
    }
private:
    unsigned long m_IterationNumber;
};

int main(int argc, const char * argv[]) {
    
    //to determine time of computation
    std::srand(time(NULL));
    std::time_t tbegin,tend;
    double texec = 0;
    tbegin = std::time(NULL);
   
    //lecture des images
    
    string filenameUS;
    string filenameIRM;
    string outputFilename;
    
    /*********************
     * RUNTIME ARG
     ******************/
    
    cout<<"runtime arguments acquisition"<<endl;
    
    for(int i = 0; i < argc; ++i)
    {
        //input image US
        if(strcmp(argv[i], "-iUS")==0)
        {
            i++;
            filenameUS = argv[i];
        }
        
        //input image IRM
        if(strcmp(argv[i], "-iIRM")==0)
        {
            i++;
            filenameIRM = argv[i];
        }
        
        if(strcmp(argv[i], "-o")==0)
        {
            i++;
            outputFilename = argv[i];
        }
        
        
    }
    
    /**********************
     * VERIFICATION INPUTS
     ********************/
    
    cout<<"input error handling"<<endl;
    
    if(filenameUS == "")
    {
        cerr<<"Input US file not provided"<<endl;
        return EXIT_FAILURE;
        
    }
    
    if(filenameIRM == "")
    {
        cerr<<"Input MRI file not provided"<<endl;
        return EXIT_FAILURE;
    }
    
    if(outputFilename == "")
    {
        cerr<<"output path not provided"<<endl;
        return EXIT_FAILURE;
    }
    
    /**********************
     * US READING
     *********************/
    
    cout<<"Reading images"<<endl;
    
    ImageType::Pointer image_US = ImageType::New();
    ImageType::Pointer image_IRM = ImageType::New();
    
    
    ReaderType::Pointer reader1 = ReaderType::New();
    itk::MetaImageIO::Pointer m_io = itk::MetaImageIO::New();
    reader1->SetImageIO(m_io);
    reader1->SetFileName(filenameUS);
    try {
        reader1->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"Error while reading US image"<<endl;
        cerr<<e<<endl;
        EXIT_FAILURE;
    }
    
    image_US = reader1->GetOutput();
    cout<<"test lecture US"<<endl;
    cout<<"dimensions US : "<<image_US->GetLargestPossibleRegion().GetSize()<<endl;
    
    /*******************
     * DOWNSAMPLING US
     ******************/
    //RESOLUTION ADAPTATION US
    ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
    shrinkFilter->SetInput(image_US);
    //porc 6 : 3,3,2/ porc 1 5,4,3
    shrinkFilter->SetShrinkFactor(0, 3);
    shrinkFilter->SetShrinkFactor(1, 3);
    shrinkFilter->SetShrinkFactor(2, 2);
    try {
        shrinkFilter->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while downsampling US image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    ImageType::Pointer US_shrunk = shrinkFilter->GetOutput();
    
    //verification ecriture de l'image
        WriterType::Pointer writer6 = WriterType::New();
        string out6 = "/Users/maximegerard/Documents/ShrunkUS.nii.gz";
      itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
        writer6->SetImageIO(io);
        writer6->SetInput(US_shrunk);
        writer6->SetFileName(out6);
        try {
            writer6->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"error while writing rescaled image"<<endl;
            cerr<<e<<endl;
            return EXIT_FAILURE;
        }
    
    cout<<"done writing shrunk US"<<endl;
    
    
    
    
//    //min max image US
//    MinMaxCalculatorType::Pointer minMaxUS = MinMaxCalculatorType::New();
//    minMaxUS->SetImage(image_US);
//    minMaxUS->Compute();
//    
//    cout<<"intensity range US image : "<<"[ "<<minMaxUS->GetMinimum()<<","<<minMaxUS->GetMaximum()<<" ]"<<endl;
    
    /**********************
     * IRM READING
     *********************/
    
    ReaderType::Pointer reader2 = ReaderType::New();
    //itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
    reader2->SetImageIO(io);
    reader2->SetFileName(filenameIRM);
    try {
        reader2->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"Error while reading IRM image"<<endl;
        cerr<<e<<endl;
        EXIT_FAILURE;
    }
    
    image_IRM = reader2->GetOutput();
    cout<<"test lecture IRM"<<endl;
    cout<<"dimensions IRM : "<<image_IRM->GetLargestPossibleRegion().GetSize()<<endl;
    
//    //min max image IRM
//    MinMaxCalculatorType::Pointer minMaxIRM = MinMaxCalculatorType::New();
//    minMaxIRM->SetImage(image_IRM);
//    minMaxIRM->Compute();
//    
//    cout<<"initial intensity range IRM image : "<<"[ "<<minMaxIRM->GetMinimum()<<","<<minMaxIRM->GetMaximum()<<" ]"<<endl;
    
    cout<<"done reading images"<<endl;
    
    /****************
    * RESCALING IRM
     ***************/
    
    cout<<"rescaling de l'image IRM"<<endl;
    RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetInput(image_IRM);
    rescaler->SetOutputMinimum(0);
    rescaler->SetOutputMaximum(255);
    try {
        rescaler->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while rescaling image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    
    ImageType::Pointer rescaled_IRM = rescaler->GetOutput();
    
//    //min max image IRM
//    MinMaxCalculatorType::Pointer minMaxIRM2 = MinMaxCalculatorType::New();
//    minMaxIRM2->SetImage(rescaled_IRM);
//    minMaxIRM2->Compute();
//    
//    cout<<"rescaled intensity range IRM image : "<<"[ "<<minMaxIRM2->GetMinimum()<<","<<minMaxIRM2->GetMaximum()<<" ]"<<endl;
    
//    WriterType::Pointer writer6 = WriterType::New();
//    string out6 = "/Users/maximegerard/Documents/IRMRescaled.nii.gz";
//    writer6->SetImageIO(io);
//    writer6->SetInput(rescaled_IRM);
//    writer6->SetFileName(out6);
//    try {
//        writer6->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error while writing rescaled image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
//   
    
    
    /*************************
     * RECALAGE
    *************************/
    
    cout<<"registration routine"<<endl;
    
    ////////////////
    //1. TRANSFORM /
    ////////////////
    
    cout<<"creation tsf"<<endl;
    
    //euler tsf
    //TransformType::Pointer transform = TransformType::New();
    //transform->SetIdentity();
    
    //tsl transform
    TranslationType::Pointer transform = TranslationType::New();
    
    /////////////////
    // 2. OPTIMIZER /
    /////////////////
    
 //LM optimizer
  //OptimizerType::Pointer optimizer = OptimizerType::New();
//    optimizer->SetUseCostFunctionGradient(false);
//    OptimizerType::ScalesType scales(transform->GetNumberOfParameters());
//    //setting the scales and ranges for the parameters of the transform
//    const double translationScale = 1000.0;
//    const double rotationScale = 1.0;
//    
//    scales[0] = 1.0;
//    scales[1] = 1.0;
//    scales[2] = 1.0;
//    scales[3] = 1.0;
//    scales[4] = 1.0;
//    scales[5] = 1.0;
//    
//    unsigned long numberOfIterations = 200;
//    double gradientTolerance = 1e-4;
//    double valueTolerance = 1e-4;
//    double epsilonFunction = 1e-5;
//    
//    optimizer->SetScales(scales);
//    optimizer->SetNumberOfIterations(numberOfIterations);
//    optimizer->SetValueTolerance(valueTolerance);
//    optimizer->SetGradientTolerance(gradientTolerance);
//    optimizer->SetEpsilonFunction(epsilonFunction);
    
    //Amoeba optimizer
    
    cout<<"creation optimizer"<<endl;
    
    OptimizerType::Pointer optimizer = OptimizerType::New();
    
    
    //interpolateur et metrique
    
    /////////////////////
    // 3. INTERPOLATOR //
    /////////////////////
    
    cout<<"creation interpolateur"<<endl;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    
    
    /////////////////
    //  4. REGISTOR /
    ////////////////
    //pipeline de recalage
    
    RegistrationType::Pointer registration = RegistrationType::New();
    
    //settings des blocs constitutifs
   
    cout<<"settings des blocs constitutifs"<<endl;
  
    registration->SetOptimizer(optimizer);
    registration->SetTransform(transform);
    registration->SetInterpolator(interpolator);

    
    ///////////////
    // 4. METRIC //
    ///////////////
    cout<<"setting metrique"<<endl;
    LC2MetricType::Pointer metric = LC2MetricType::New();
    
    registration->SetMetric(metric);
    
    
    //setting des images
    cout<<"setting images"<<endl;
    registration->SetFixedImage(US_shrunk);
    registration->SetMovingImage(rescaled_IRM);
    //just to be sure for the metric
    metric->SetFixed(US_shrunk);
    metric->ComputeMask();
    metric->SetMoving(rescaled_IRM);
    registration->SetFixedImageRegion(US_shrunk->GetBufferedRegion());
    
    //setting de la tsf initiale à 0,0,0 de tsl
    transform->SetIdentity();
    
    RegistrationType::ParametersType initialParameters = transform->GetParameters();
    initialParameters[0] = 0.0;
    initialParameters[1] = 0.0;
    initialParameters[2] = 0.0;
    
    registration->SetInitialTransformParameters(initialParameters);
    
    cout<<"Initial transform param = "<<initialParameters<<endl;
    
    //setting des params de l'optimizer
    const unsigned int numberOfParameters = transform->GetNumberOfParameters();
    OptimizerType::ParametersType simplexDelta(numberOfParameters);
    simplexDelta.fill(5.0);
    
    optimizer->AutomaticInitialSimplexOff();
    optimizer->SetInitialSimplexDelta(simplexDelta);
    
    //tolerance setting on param units and LC2 units
    optimizer->SetParametersConvergenceTolerance(0.1);
    optimizer->SetFunctionConvergenceTolerance(0.001);
    
    //max number of iterations to stop if no convergence
    optimizer->SetMaximumNumberOfIterations(200);
    
    //on le met en maximisation
    optimizer->SetMaximize(true);
    
    //Command observer
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    optimizer->AddObserver(itk::IterationEvent(), observer);
    

    cout<<"registration"<<endl;
    
        try {
            registration->Update();
            cout<<"Optimizer stop condition : "<<registration->GetOptimizer()->GetStopConditionDescription()<<endl;
        } catch (itk::ExceptionObject &e) {
            cerr<<"erreur dans le recalage"<<endl;
            cerr<<e<<endl;
            return EXIT_FAILURE;
        }
    

    RegistrationType::ParametersType finalParameters = registration->GetLastTransformParameters();
    const double TX = finalParameters[0];
    const double TY = finalParameters[1];
    const double TZ = finalParameters[2];
    
    double bestValue = optimizer->GetValue();
    
    //Print out results
    
    cout<<"Results : "<<endl;
    cout<<"Translation vector : "<<"[ "<<TX<<", "<<TY<<", "<<TZ<<" ]"<<endl;
    cout<<"metric value : "<<bestValue<<endl;
    
    cout<<"formation de l'image recalee"<<endl;
    
    TranslationType::Pointer finalTsl = TranslationType::New();
    finalTsl->SetParameters(finalParameters);
    finalTsl->SetFixedParameters(transform->GetFixedParameters());
    
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetTransform(finalTsl);
    resampler->SetInput(rescaled_IRM);
    resampler->SetOutputDirection(rescaled_IRM->GetDirection());
    resampler->SetSize(rescaled_IRM->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputSpacing(rescaled_IRM->GetSpacing());
    ImageType::PointType origine = rescaled_IRM->GetOrigin();
    origine[0] = origine[0] - TX;
    origine[1] = origine[1] - TY;
    origine[2] = origine[2] - TZ;
    resampler->SetOutputOrigin(origine);
    resampler->SetDefaultPixelValue(0);
    
    ImageType::Pointer outputImage = ImageType::New();
    outputImage = resampler->GetOutput();
    
    WriterType::Pointer writer7 = WriterType::New();
    string out7 ="/Users/maximegerard/Documents/finalregistreredIRM.nii.gz";
    writer7->SetImageIO(io);
    writer7->SetInput(outputImage);
    writer7->SetFileName(out7);
    try {
        writer7->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error whilte writing registered image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    
    
    
    
    //TEST METRIC LC2
   // cout<<"test metric"<<endl;
//    metric->SetFixed(US_shrunk);
//    metric->ComputeMask();
////    metric->SetMoving(rescaled_IRM);
////    metric->test();
//    
//
//    
//    
//    //TEST FORCE BRUTE
//    cout<<"test force brute"<<endl;
//    ResampleFilterType::Pointer resamplefilter = ResampleFilterType::New();
//    resamplefilter->SetInput(rescaled_IRM);
//    resamplefilter->SetOutputSpacing(rescaled_IRM->GetSpacing());
//    resamplefilter->SetSize(rescaled_IRM->GetLargestPossibleRegion().GetSize());
//    resamplefilter->SetOutputDirection(rescaled_IRM->GetDirection());
//    ImageType::PointType origine = rescaled_IRM->GetOrigin();
//    ImageType::Pointer transformed_IRM;
//    ImageType::Pointer final_MRI = ImageType::New();
//    cout<<"origine before"<<origine<<endl;
//    double max =0;
//    
//    cout<<"test verification"<<endl;
//    TranslationType::OutputVectorType best_tsl;
//    
//
//    for(double i=-7;i<8;i++)
//    {
//        TranslationType::Pointer tsf = TranslationType::New();
//        TranslationType::OutputVectorType tsl;
// 
//        
//        tsl[0]=i;
//        tsl[1]=0;
//        tsl[2]=0;
//        cout<<"tsl value : "<<tsl<<endl;
//        
//        tsf->Translate(tsl);
//        resamplefilter->SetTransform(tsf);
//        origine = rescaled_IRM->GetOrigin();
//        cout<<"origine before : "<<origine<<endl;
//        origine[0] = origine[0]-tsl[0];
//        origine[1] = origine[1]-tsl[1];
//        origine[2] = origine[2]-tsl[2];
//        cout<<"origine after: "<<origine<<endl;
//        resamplefilter->SetOutputOrigin(origine);
//        try {
//            resamplefilter->Update();
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"error while translating image"<<endl;
//            cerr<<e<<endl;
//            return EXIT_FAILURE;
//        }
//        
//        ImageType::Pointer result = resamplefilter->GetOutput();
//        metric->SetMoving(result);
//        double score = metric->test();
//        if(score>max)
//        {
//            max=score;
//            best_tsl=tsl;
//            
//        }
//        
//    }
//    
//    cout<<"best tsl : "<<best_tsl<<endl;
//    
//    TranslationType::Pointer final_tsf = TranslationType::New();
//    final_tsf->Translate(best_tsl);
//    resamplefilter->SetTransform(final_tsf);
//    origine = rescaled_IRM->GetOrigin() - best_tsl;
//    resamplefilter->SetOutputOrigin(origine);
//    try {
//        resamplefilter->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error while translating image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
//    
//  final_MRI = resamplefilter->GetOutput();
//    
//    
//
//
//    
//    
//        //write image result
//        WriterType::Pointer writerfinal = WriterType::New();
//        string outfinal ="/Users/maximegerard/Documents/finalIRM.nii.gz"  ;
//        writerfinal->SetImageIO(io);
//        writerfinal->SetFileName(outfinal);
//        writerfinal->SetInput(final_MRI);
//        try {
//            writerfinal->Update();
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"error whilte writing registered image"<<endl;
//            cerr<<e<<endl;
//            return EXIT_FAILURE;
//        }
//    
    cout<<"done writing translation image"<<endl;
    
    
    
    
    
    

    
    std::cout << "Done running routine!\n";
  
    tend = std::time(NULL);
    texec = std::difftime(tend,tbegin);
    
    cout<<"temps en s : "<<texec<<endl;
    return 0;
}
