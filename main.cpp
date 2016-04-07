//
//  main.cpp
//  LC2_M
// Calcul de la metrique LC2
//  Created by Maxime Gérard on 17/12/15.
//  Copyright © 2015 Maxime Gérard. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>

#include "gradient.h"
#include "LC2_function.hpp"



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
#include "itkAffineTransform.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkCastImageFilter.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkLC2MVImageToImageMetric.h"
#include "itkBsplineTransform.h"
#include "itkAdaptiveHistogramEqualizationImageFilter.h"

#include <cmath>

#include <dlib/optimization.h>
#include <dlib/matrix.h>

using namespace std;
using namespace dlib;




//images
typedef itk::Image<double,3> ImageType;
typedef itk::Image<unsigned char, 3> BinaryImageType;
typedef itk::Image<double,2> Image2DType;

//IO
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageFileReader<Image2DType> Reader2DType;
typedef itk::ImageFileReader<BinaryImageType> BinaryReaderType;
typedef itk::ImageFileWriter<ImageType> WriterType;

//iteration over images
typedef itk::ImageRegionIterator<ImageType> IteratorType;

//recalage
typedef itk::AmoebaOptimizer AmoebaOptimizerType;
typedef itk::LevenbergMarquardtOptimizer LMOptimizerType;
typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;
typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
typedef itk::Euler3DTransform<double> EulerTransformType;
typedef itk::LC2ImageToImageMetric<ImageType, ImageType> LC2MetricType;
typedef itk::TranslationTransform<double,3> TranslationType;
typedef itk::AffineTransform<double,3> AffineTransformType;
typedef itk::BSplineTransform<double,3,3> BSplineTransformType;
typedef itk::AdaptiveHistogramEqualizationImageFilter<ImageType> HistoEqualizerType;

typedef itk::LC2MVImageToImageMetric<ImageType, ImageType> LC2MVMetricType;

//pour le mapping de l'image mobile registree ac tsf determinee par registration framework
typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;

//pour la mise a l'echelle de la bande passante des intensites
typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> RescaleFilterType;
typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxCalculatorType;
typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkFilterType;

//MultiResolutionRegistration


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

////classe pour suivre le recalage multires
//template <typename TRegistration>
//class RegistrationInterfaceCommand : public itk::Command
//{
//public:
//    typedef  RegistrationInterfaceCommand   Self;
//    typedef  itk::Command                   Superclass;
//    typedef  itk::SmartPointer<Self>        Pointer;
//    itkNewMacro( Self );
//protected:
//    RegistrationInterfaceCommand() {};
//public:
//    typedef   TRegistration                              RegistrationType;
//    typedef   RegistrationType *                         RegistrationPointer;
//    typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType;
//    typedef   OptimizerType *                            OptimizerPointer;
//    void Execute(itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
//    {
//        if( !(itk::IterationEvent().CheckEvent( &event )) )
//        {
//            return;
//        }
//        RegistrationPointer registration = static_cast<RegistrationPointer>( object );
//        if(registration == ITK_NULLPTR)
//        {
//            return;
//        }
//        OptimizerPointer optimizer = static_cast< OptimizerPointer >(registration->GetModifiableOptimizer() );
//        std::cout << "-------------------------------------" << std::endl;
//        std::cout << "MultiResolution Level : "
//        << registration->GetCurrentLevel()  << std::endl;
//        std::cout << std::endl;
//        if ( registration->GetCurrentLevel() == 0 )
//        {
//            optimizer->SetMaximumStepLength( 16.00 );
//            optimizer->SetMinimumStepLength( 0.01 );
//        }
//        else
//        {
//            optimizer->SetMaximumStepLength( optimizer->GetMaximumStepLength() / 4.0 );
//            optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength() / 10.0 );
//        }
//    }
//    void Execute(const itk::Object * , const itk::EventObject & ) ITK_OVERRIDE
//    { return; }
//};
//
//class CommandIterationUpdate : public itk::Command
//{
//public:
//    typedef  CommandIterationUpdate   Self;
//    typedef  itk::Command             Superclass;
//    typedef  itk::SmartPointer<Self>  Pointer;
//    itkNewMacro( Self );
//protected:
//    CommandIterationUpdate() {};
//public:
//    typedef   itk::RegularStepGradientDescentOptimizer OptimizerType;
//    typedef   const OptimizerType *                    OptimizerPointer;
//    void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
//    {
//        Execute( (const itk::Object *)caller, event);
//    }
//    void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
//    {
//        OptimizerPointer optimizer = static_cast< OptimizerPointer >( object );
//        if( !(itk::IterationEvent().CheckEvent( &event )) )
//        {
//            return;
//        }
//        std::cout << optimizer->GetCurrentIteration() << "   ";
//        std::cout << optimizer->GetValue() << "   ";
//        std::cout << optimizer->GetCurrentPosition() << std::endl;
//    }
//};

int main(int argc, const char * argv[]) {
    
    //to determine time of computation
    std::srand(time(NULL));
    std::time_t tbegin,tend;
    double texec = 0;
    tbegin = std::time(NULL);
   
    //lecture des images
    
    string filenameUS;
    string filenameIRM;
    string outputPath;
    string filenameMaskLiver;
    
    bool useLiverMask=false;
    
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
        
        if(strcmp(argv[i], "-iMaskLiver")==0)
        {
            i++;
            filenameMaskLiver= argv[i];
            useLiverMask=true;
            cout<<"Use of mask image of liver"<<endl;
        }
        
        if(strcmp(argv[i], "-o")==0)
        {
            i++;
            outputPath = argv[i];
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
    
    
    if(outputPath == "")
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
    
    
    //min max image US
    MinMaxCalculatorType::Pointer minMaxUS = MinMaxCalculatorType::New();
    minMaxUS->SetImage(image_US);
    minMaxUS->Compute();
    
    cout<<"intensity range US image : "<<"[ "<<minMaxUS->GetMinimum()<<","<<minMaxUS->GetMaximum()<<" ]"<<endl;
    
    /**********************
     * IRM READING
     *********************/
    
    ReaderType::Pointer reader2 = ReaderType::New();
    itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
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
    
    //min max image IRM
    MinMaxCalculatorType::Pointer minMaxIRM = MinMaxCalculatorType::New();
    minMaxIRM->SetImage(image_IRM);
    minMaxIRM->Compute();
    
    cout<<"initial intensity range IRM image : "<<"[ "<<minMaxIRM->GetMinimum()<<","<<minMaxIRM->GetMaximum()<<" ]"<<endl;
    
    cout<<"done reading images"<<endl;
    
    /********************
     * Liver mask reading
     **********************/
    
    //gotta declare it outside if loop to be visible within the entire main function
    BinaryImageType::Pointer LiverMask = BinaryImageType::New();
    
    if(useLiverMask)
    {
        cout<<"Reading MRI liver mask"<<endl;
        BinaryReaderType::Pointer BinReader = BinaryReaderType::New();
        BinReader->SetImageIO(io);
        BinReader->SetFileName(filenameMaskLiver);
        try {
            BinReader->Update();
        } catch (itk::ExceptionObject &e) {
            cerr<<"Error while reading liver mask"<<endl;
            cerr<<e<<endl;
        }
        
       LiverMask =BinReader->GetOutput();
    }

    
    
    /*******************
     * DOWNSAMPLING US
     ******************/
    //RESOLUTION ADAPTATION US
    ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
    shrinkFilter->SetInput(image_US);
    
    //recuperation des spacing des deux images pour savoir quel facteur utiliser
    ImageType::SpacingType spacingUS = image_US->GetSpacing();
    ImageType::SpacingType spacingIRM = image_IRM->GetSpacing();
    //calcul du spacing factor
    int shrinkX = int(2*spacingIRM[0]/spacingUS[0]);
    int shrinkY = int(2*spacingIRM[1]/spacingUS[1]);
    int shrinkZ = int(2*spacingIRM[2]/spacingUS[2]);
    
    cout<<"shrinking factors : "<<shrinkX<<", "<<shrinkY<<", "<<shrinkZ<<endl;

    //porc 6 : (3,3,2)/ porc 1 (5,4,3)
    shrinkFilter->SetShrinkFactor(0, shrinkX);
    shrinkFilter->SetShrinkFactor(1, shrinkY);
    shrinkFilter->SetShrinkFactor(2, shrinkZ);
    try {
        shrinkFilter->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while downsampling US image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
   ImageType::Pointer US_shrunk = shrinkFilter->GetOutput();

//        //verification ecriture de l'image
//            WriterType::Pointer writer6 = WriterType::New();
//            string out6 = outputPath+"/ShrunkUS.nii.gz";
//            //itk::NiftiImageIO::Pointer io = itk::NiftiImageIO::New();
//            writer6->SetImageIO(io);
//            writer6->SetInput(US_shrunk);
//            writer6->SetFileName(out6);
//            try {
//                writer6->Update();
//            } catch (itk::ExceptionObject &e) {
//                cerr<<"error while writing rescaled image"<<endl;
//                cerr<<e<<endl;
//                return EXIT_FAILURE;
//            }
//
    cout<<"done writing shrunk US"<<std::endl;
    cout<<"taille US shrunk : "<<US_shrunk->GetLargestPossibleRegion().GetSize()<<std::endl;
    
//    /**********************
//     * RESCALING US
//     **************************/
//    //put US image in the MRI intensity bandwidth
//    
//    cout<<"Rescaling US image intensities"<<endl;
//    
//    RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
//    rescaler->SetInput(US_shrunk);
//    //rescaler->SetInput(MRI_shrunk);
//    rescaler->SetOutputMinimum(minMaxIRM->GetMinimum());
//    rescaler->SetOutputMaximum(minMaxIRM->GetMaximum());
//    try {
//        rescaler->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error while rescaling image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
//    
//    
//    ImageType::Pointer rescaled_US= rescaler->GetOutput();
//    
//    //writer result
//    WriterType::Pointer writer11 = WriterType::New();
//    string out6 = outputPath+"/rescaled_US.nii.gz";
//    writer11->SetImageIO(io);
//    writer11->SetInput(rescaled_US);
//    writer11->SetFileName(out6);
//    try {
//        writer11->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error while writing rescaled image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
    
    /****************
     * RESCALING IRM
     ***************/
    
    //MAKE SUR MRI INTENSITIES ARE WITHIN [0,255] -> why not put the US data in the MRI range ?
    cout<<"rescaling de l'image IRM"<<endl;
    RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetInput(image_IRM);
    //rescaler->SetInput(MRI_shrunk);
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
    
    //ecrire res
    WriterType::Pointer writer10 = WriterType::New();
    writer10->SetImageIO(io);
    writer10->SetInput(rescaled_IRM);
    string out10 = outputPath+"/rescaled_MRI.nii.gz";
    writer10->SetFileName(out10);
    try {
        writer10->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"Error while writing rescaled MRI"<<endl;
        cerr<<e<<endl;
    }

    //min max image IRM
    MinMaxCalculatorType::Pointer minMaxIRM2 = MinMaxCalculatorType::New();
    minMaxIRM2->SetImage(rescaled_IRM);
    minMaxIRM2->Compute();

    cout<<"rescaled intensity range IRM image : "<<"[ "<<minMaxIRM2->GetMinimum()<<","<<minMaxIRM2->GetMaximum()<<" ]"<<endl;

    
    /*******************
     * DOWNSAMPLING IRM
     *******************/
    
    //APPROCHE MUTLIRES
    
    ShrinkFilterType::Pointer shrinkerMRI = ShrinkFilterType::New();
    shrinkerMRI->SetInput(image_IRM);
    shrinkerMRI->SetShrinkFactor(0, 2); //2 pour MRI
    shrinkerMRI->SetShrinkFactor(1, 2);
    shrinkerMRI->SetShrinkFactor(2, 2);
    try {
        shrinkerMRI->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error while shrinking MRI image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    ImageType::Pointer MRI_shrunk = shrinkerMRI->GetOutput();
    
    cout<<"taille MRI shrunk : "<<MRI_shrunk->GetLargestPossibleRegion().GetSize()<<endl;
    
//    /********************************
//     * Egalisation d'hisogramme IRM
//     *******************************/
//    
//    cout<<"Histogramme equalization MRI image"<<endl;
//    
//    HistoEqualizerType::Pointer equalizer = HistoEqualizerType::New();
//    equalizer->SetInput(rescaled_IRM);
//    //equalizer->SetRadius(1);
//    equalizer->SetAlpha(0);
//    equalizer->SetBeta(0.5);
//    try {
//        equalizer->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"Error while enhancing IRM image contrast"<<endl;
//        cerr<<e<<endl;
//    }
//    
//    ImageType::Pointer IRM_enhanced = equalizer->GetOutput();
//    
//    WriterType::Pointer writer9 = WriterType::New();
//    writer9->SetImageIO(io);
//    writer9->SetInput(IRM_enhanced);
//    string out9 = outputPath+"/equalizedMRI.nii.gz";
//    writer9->SetFileName(out9);
//    try {
//        writer9->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"Error while writing equalized MRI"<<endl;
//        cerr<<e<<endl;
//    }
    
    

    
    
    /*************************
     * RECALAGE
    *************************/
    
    cout<<"registration routine"<<endl;
    
    //test cropping MRI
    
//    LC2_function LC2 = LC2_function(rescaled_IRM, US_shrunk);
//    LC2.limitMRI();
    
    //they must be column_vector !!!
    
    //BOBYQA initiaux FULL RES
    //Porc 1
    //best score : 0.183876
    //final parameters : [-0.03514169419515513, 0.04780031905516767, 0.08510374453376125, -3.16002229212777, -0.7603111683282886, -1.5771778577960676]
    
    //porc 2
    //MRI Shrunk : best score 1ere etape : 0.149085
    //final parameters : [-0.006769331840125886, 0.09820386001984022, -0.01826162744483251, 2.7224875212753967, 0.7292929664137704, 3.0486718768199808]
    
    //full resolution MRI
//    best score 1ere etape : 0.156658
//    final parameters : [0.0036297748801233246, 0.013143829260491974, 0.0010727709811333918, 0.08008875517288504, -0.04607054405119398, 1.0090378953697132]
    
    //MRI rescaled
    
    //porc4
//    best score : 0.185428
//    final parameters : [-0.025424753944333096, 0.39962575090364716, 0.013544113893079517, -0.047444368532859084, -0.739287807247723, -0.5011674690804339]
    
    //porc6
//    best score : 0.150082
//    final parameters : [-0.0150802635531933, -0.010100620088642457, 0.004725328308930741, 8.30720384677868, 0.3776652489827957, 0.9237480980253407]


    //MRI Shrunk + mask
//    best score 1ere etape : 0.195807
//    final parameters : [-0.004576003822187878, -0.018662112915280166, 0.021674120903660443, 1.335316780176604, 0.012001995175122276, -0.08281897027573802]
    
    /**********
     *IDEAS TO TEST
     ************/

    //using mask of liver to limit ROI
    //try with arterial phased MRA (p6)
    //try playing on parameters to understand BOBYQA better, what about CMA-ES ?
    //GPU programming ?
  
    
    matrix<double> initialParameters (6,1);
    initialParameters(0) = 0;
    initialParameters(1) = 0;
    initialParameters(2) = 0;
    initialParameters(3) = 0;
    initialParameters(4) = 0;
    initialParameters(5) = 0;
    
    
    //cout<<"test params intial : "<<initialParameters<<endl;
    
    matrix<double> x_lower (6,1); //-0.5 rot, -10 trans
    x_lower(0) = -1;
    x_lower(1) = -1;
    x_lower(2) = -1;
    x_lower(3) = -1;
    x_lower(4) = -1;
    x_lower(5) = -1;
    
    matrix<double> x_upper (6,1); //0.5 rot, 10 trans
    x_upper(0) = 1;
    x_upper(1) = 1;
    x_upper(2) = 1;
    x_upper(3) = 1;
    x_upper(4) = 1;
    x_upper(5) = 1;
    
    double m = 13;
    
    //rho begin
    double radius = 0.9; //0.9 //zero makes no sense !!!
    
    //rho end
    double precision = 0.1; //0.01 0.001
    
    //niter
    double nombreIteration = 200;
    
    LC2_function LC2 = LC2_function(MRI_shrunk,US_shrunk);//MRI_shrunk or rescaled
    //make sure that the mask is computed on US_Shrunk but that we use the rescaled image to compute LC2
    //LC2.setMovingImage(rescaled_US);
    LC2.setMaxRot(0.3);
    LC2.setMaxTrans(10);
    LC2.setRadius(radius);
    LC2.setOutputPath(outputPath);
    
    if(useLiverMask) LC2.setLiverMask(LiverMask);
    
   double best_score = find_max_bobyqa(LC2, initialParameters, m, x_lower, x_upper, LC2.getRadius(), precision, nombreIteration);
    //double best_score = find_min_bobyqa(LC2_function(rescaled_IRM, US_shrunk), initialParameters, m, x_lower, x_upper, radius, precision, nombreIteration);
    
    cout<<"best score 1ere etape : "<<best_score<<endl;
    //cout<<"rough parameters : "<<initialParameters<<endl;

//    EulerTransformType::ParametersType Step1Parameters(6);
//    
//    Step1Parameters[0] = initialParameters(0)*LC2.getMaxRot()/(LC2.getRadius());
//    Step1Parameters[1]= initialParameters(1)*LC2.getMaxRot()/(LC2.getRadius());
//    Step1Parameters[2] = initialParameters(2)*LC2.getMaxRot()/(LC2.getRadius());
//    Step1Parameters[3] = initialParameters(3)*LC2.getMaxTrans()/(LC2.getRadius());
//    Step1Parameters[4] = initialParameters(4)*LC2.getMaxTrans()/(LC2.getRadius());
//    Step1Parameters[5] = initialParameters(5)*LC2.getMaxTrans()/(LC2.getRadius());
//    cout<<"test best parameters 1ere etape: "<<Step1Parameters<<endl;
    
    
    //enregistrement resultats
    EulerTransformType::ParametersType finalParameters(6);
    finalParameters[0] = initialParameters(0)*LC2.getMaxRot()/(LC2.getRadius());
    finalParameters[1] = initialParameters(1)*LC2.getMaxRot()/(LC2.getRadius());
    finalParameters[2] = initialParameters(2)*LC2.getMaxRot()/(LC2.getRadius());
    finalParameters[3] = initialParameters(3)*LC2.getMaxTrans()/(LC2.getRadius());
    finalParameters[4] = initialParameters(4)*LC2.getMaxTrans()/(LC2.getRadius());
    finalParameters[5] = initialParameters(5)*LC2.getMaxTrans()/(LC2.getRadius());
    
//    //transformation de l'image US pour la reutiliser comme nouveau point de départ avec une plus grande resolution IRM
//    
//    EulerTransformType::Pointer Step1Tsf = EulerTransformType::New();
//    Step1Tsf->SetParameters(Step1Parameters);
//    
//    cout<<"step 1 parameters : "<<Step1Tsf->GetParameters()<<endl;
//    
//    ImageType::SizeType sizeUS = US_shrunk->GetLargestPossibleRegion().GetSize();
//    ImageType::PointType origin = US_shrunk->GetOrigin();
//    ImageType::SpacingType spacing = US_shrunk->GetSpacing();
//    ImageType::PointType center;
//    center[0] = origin[0]+spacing[0]*sizeUS[0]/2;
//    center[1] = origin[1]+spacing[1]*sizeUS[1]/2;
//    center[2] = origin[2]+spacing[2]*sizeUS[2]/2;
//    
//    
//    EulerTransformType::ParametersType eulerFixedParameters(3);
//    eulerFixedParameters[0] =center[0];
//    eulerFixedParameters[1] =center[1];
//    eulerFixedParameters[2] =center[2];
//    
//   Step1Tsf->SetFixedParameters(eulerFixedParameters);
//    
//    //ici c'est toujours le US_shrunk qu'on utilise puisqu'on est pas encore a l'etape finale
//    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
//    resampler->SetTransform(Step1Tsf);
//    resampler->SetInput(US_shrunk);
//    resampler->SetOutputDirection(Step1Tsf->GetInverseMatrix()*US_shrunk->GetDirection());
//    resampler->SetSize(US_shrunk->GetLargestPossibleRegion().GetSize());
//    resampler->SetOutputSpacing(US_shrunk->GetSpacing());
//    resampler->SetOutputOrigin(Step1Tsf->GetInverseTransform()->TransformPoint(US_shrunk->GetOrigin()));
//    
//    ImageType::Pointer US_step1 = ImageType::New();
//    US_step1 = resampler->GetOutput(); //US_step1 est notre nouvelle image de travail
//    //apres on peut faire une comporsition de transformmee
//    
//    //on applique aussi a l'image de base qui est celle qu'on veut au final tsf
//    resampler->SetInput(image_US);
//    resampler->SetOutputDirection(Step1Tsf->GetInverseMatrix()*image_US->GetDirection());
//    resampler->SetSize(image_US->GetLargestPossibleRegion().GetSize());
//    resampler->SetOutputSpacing(image_US->GetSpacing());
//    resampler->SetOutputOrigin(Step1Tsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin()));
//    
//    ImageType::Pointer temp_image_US = ImageType::New();
//    temp_image_US = resampler->GetOutput();
//    
//    
//    WriterType::Pointer writer7 = WriterType::New();
//    string out7 =outputPath"/IntermediaryregistreredUSBOBYQA.nii.gz";
//    writer7->SetImageIO(io);
//    writer7->SetInput(temp_image_US);
//    writer7->SetFileName(out7);
//    try {
//        writer7->Update();
//    } catch (itk::ExceptionObject &e) {
//        cerr<<"error whilte writing registered image"<<endl;
//        cerr<<e<<endl;
//        return EXIT_FAILURE;
//    }
    
    
//   
//    /*****************************
//     *BOBYQA 2
//     *************************/
//    cout<<"Nouvelle etape BOBYQA"<<endl;
//    //mtn on utiliser IRM rescaled et US_step1
//    
//    matrix<double> initialParameters2 (6,1);
//    initialParameters2(0) = 0;
//    initialParameters2(1) = 0;
//    initialParameters2(2) = 0;
//    initialParameters2(3) = 0;
//    initialParameters2(4) = 0;
//    initialParameters2(5) = 0;
//    
//    
//    //cout<<"test params intial : "<<initialParameters<<endl;
//    
//    matrix<double> x_lower2 (6,1); //-0.5 rot, -10 trans
//    x_lower2(0) = -1;
//    x_lower2(1) = -1;
//    x_lower2(2) = -1;
//    x_lower2(3) = -1;
//    x_lower2(4) = -1;
//    x_lower2(5) = -1;
//    
//    matrix<double> x_upper2 (6,1); //0.5 rot, 10 trans
//    x_upper2(0) = 1;
//    x_upper2(1) = 1;
//    x_upper2(2) = 1;
//    x_upper2(3) = 1;
//    x_upper2(4) = 1;
//    x_upper2(5) = 1;
//    
//    double m2 = 13;
//    
//    //rho begin
//    double radius2 = 0.5; //zero makes no sense !!!
//    
//    //rho end
//    double precision2 = 0.001;
//    
//    //niter
//    double nombreIteration2 = 200;
//    
//    LC2_function LC2_2 = LC2_function(MRI_shrunk, US_step1);
//    LC2_2.setMaxRot(0.25);
//    LC2_2.setMaxTrans(5);
//    LC2_2.setRadius(radius2);
//    
//    double best_score2 = find_max_bobyqa(LC2_2, initialParameters2, m2, x_lower2, x_upper2, LC2_2.getRadius(), precision2, nombreIteration2);
//    //double best_score = find_min_bobyqa(LC2_function(rescaled_IRM, US_shrunk), initialParameters, m, x_lower, x_upper, radius, precision, nombreIteration);
//    
//    cout<<"best score 2eme etape : "<<best_score2<<endl;
//  
//    
//    EulerTransformType::ParametersType finalParameters(6);
//    finalParameters[0] = initialParameters2(0)*LC2_2.getMaxRot()/(LC2_2.getRadius());
//    finalParameters[1] = initialParameters2(1)*LC2_2.getMaxRot()/(LC2_2.getRadius());
//    finalParameters[2] = initialParameters2(2)*LC2_2.getMaxRot()/(LC2_2.getRadius());
//    finalParameters[3] = initialParameters2(3)*LC2_2.getMaxTrans()/(LC2_2.getRadius());
//    finalParameters[4] = initialParameters2(4)*LC2_2.getMaxTrans()/(LC2_2.getRadius());
//    finalParameters[5] = initialParameters2(5)*LC2_2.getMaxTrans()/(LC2_2.getRadius());
//    
//      cout<<"best parameters 2 : "<<finalParameters<<endl;
//    
//    EulerTransformType::Pointer Step2Tsf = EulerTransformType::New();
//    Step2Tsf->SetParameters(finalParameters);
    

    
//    ImageType::SizeType sizeUS2 = US_step1->GetLargestPossibleRegion().GetSize();
//    ImageType::PointType origin2 = US_step1->GetOrigin();
//    ImageType::SpacingType spacing2 = US_step1->GetSpacing();
//    ImageType::PointType center2;
//    center2[0] = origin2[0]+spacing2[0]*sizeUS2[0]/2;
//    center2[1] = origin2[1]+spacing2[1]*sizeUS2[1]/2;
//    center2[2] = origin2[2]+spacing2[2]*sizeUS2[2]/2;
    

//    
//    EulerTransformType::ParametersType eulerFixedParameters2(3);
//    eulerFixedParameters2[0] =center2[0];
//    eulerFixedParameters2[1] =center2[1];
//    eulerFixedParameters2[2] =center2[2];
//    
//    Step2Tsf->SetFixedParameters(eulerFixedParameters2);
    
    //FINAL TSF
    
    EulerTransformType::Pointer finalTsf = EulerTransformType::New();
    finalTsf->SetParameters(finalParameters);
    
    cout<<"final parameters : "<<finalTsf->GetParameters()<<endl;
    
    //ecritire des parametres dans un fichier txt
    
    string outP =outputPath+"/parameters.txt";
    ofstream fichier(outP.c_str(),ios::out | ios::trunc );
    
    if(fichier)
    {
        fichier<<"Parameters for rigid transform : "<<endl;
        fichier<<finalTsf->GetParameters()<<endl;
        fichier<<" Score for this position : "<<endl;
        fichier<<best_score<<endl;
        fichier.close();
    }
    
    else
    {
        cerr<<"Error in opening txt file for parameters"<<endl;
    }
    
    cout<<"Writing final registered US image"<<endl;
    
    ImageType::SizeType sizeUS2 = US_shrunk->GetLargestPossibleRegion().GetSize();
    ImageType::PointType origin2 = US_shrunk->GetOrigin();
    ImageType::SpacingType spacing2 = US_shrunk->GetSpacing();
    ImageType::PointType center2;
    center2[0] = origin2[0]+spacing2[0]*sizeUS2[0]/2;
    center2[1] = origin2[1]+spacing2[1]*sizeUS2[1]/2;
    center2[2] = origin2[2]+spacing2[2]*sizeUS2[2]/2;
    
    
    EulerTransformType::ParametersType eulerFixedParameters2(3);
    eulerFixedParameters2[0] =center2[0];
    eulerFixedParameters2[1] =center2[1];
    eulerFixedParameters2[2] =center2[2];
    
    finalTsf->SetFixedParameters(eulerFixedParameters2);
    
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput(image_US);
    resampler->SetTransform(finalTsf);
    resampler->SetSize(image_US->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputSpacing(image_US->GetSpacing());
    resampler->SetOutputDirection(finalTsf->GetInverseMatrix()*image_US->GetDirection());
    resampler->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin()));
    
    ImageType::Pointer finalImage = ImageType::New();
    finalImage = resampler->GetOutput();
    
    cout<<"writing final result"<<endl;
    
    WriterType::Pointer writer8 = WriterType::New();
    string out8 =outputPath+"/finalregistreredUSBOBYQA.nii.gz";
    writer8->SetImageIO(io);
    writer8->SetInput(finalImage);
    writer8->SetFileName(out8);
    try {
        writer8->Update();
    } catch (itk::ExceptionObject &e) {
        cerr<<"error whilte writing registered image"<<endl;
        cerr<<e<<endl;
        return EXIT_FAILURE;
    }
    
    
    
    
    

    
    
//    
//    ////////////////
//    //1. TRANSFORM /
//    ////////////////
//    
//    cout<<"creation tsf"<<endl;
//    
//    //euler tsf
//    EulerTransformType::Pointer transform = EulerTransformType::New();
//    transform->SetIdentity();
//    
//    //tsl transform
//    //TranslationType::Pointer transform = TranslationType::New();
//    
//    
//    /////////////////
//    // 2. OPTIMIZER /
//    /////////////////
//    
// //LM optimizer
////    cout<<"creation LM optimizer"<<endl;
////    LMOptimizerType::Pointer optimizer = LMOptimizerType::New();
////    optimizer->SetUseCostFunctionGradient(false);
//    
//    //Amoeba optimizer
//    
//    cout<<"creation optimizer"<<endl;
//    
//    AmoebaOptimizerType::Pointer optimizer = AmoebaOptimizerType::New();
//    
//    
//    //interpolateur et metrique
//    
//    /////////////////////
//    // 3. INTERPOLATOR //
//    /////////////////////
//    
//    cout<<"creation interpolateur"<<endl;
//    InterpolatorType::Pointer interpolator = InterpolatorType::New();
//    
//    
//    /////////////////
//    //  4. REGISTOR /
//    ////////////////
//    //pipeline de recalage
//    
//    RegistrationType::Pointer registration = RegistrationType::New();
//    
//    //settings des blocs constitutifs
//   
//    cout<<"settings des blocs constitutifs"<<endl;
//  
//    registration->SetOptimizer(optimizer);
//    registration->SetTransform(transform);
//    registration->SetInterpolator(interpolator);
//
////
////    ///////////////
////    // 4. METRIC //
////    ///////////////
//    cout<<"setting metrique"<<endl;
//    LC2MetricType::Pointer metric = LC2MetricType::New();
//    
////    //test metric mutliple valued
////    LC2MVMetricType::Pointer metricMV = LC2MVMetricType::New();
////    metricMV->SetMoving(US_shrunk);
////    metricMV->SetFixed(rescaled_IRM);
////    metricMV->ComputeGradImage();
////    metricMV->ComputeMask();
////    
////    EulerTransformType::ParametersType parameters(6);
////    parameters[0]=0;
////    parameters[1]=0;
////    parameters[2]=0;
////    parameters[3]=0;
////    parameters[4]=0;
////    parameters[5]=0;
////    transform->SetParameters(parameters);
////    
////    //fixed param
////    ImageType::SizeType sizeUS = US_shrunk->GetLargestPossibleRegion().GetSize();
////    ImageType::PointType origin = US_shrunk->GetOrigin();
////    ImageType::SpacingType spacing = US_shrunk->GetSpacing();
////    ImageType::PointType center;
////    center[0] = origin[0]+spacing[0]*sizeUS[0]/2;
////    center[1] = origin[1]+spacing[1]*sizeUS[1]/2;
////    center[2] = origin[2]+spacing[2]*sizeUS[2]/2;
////    
////    
////    EulerTransformType::ParametersType eulerFixedParameters(3);
////    eulerFixedParameters[0] =center[0];
////    eulerFixedParameters[1] =center[1];
////    eulerFixedParameters[2] =center[2];
////    
////    transform->SetFixedParameters(eulerFixedParameters);
////    
////    metricMV->SetTransform(transform);
////    cout<<"params initiaux : "<<metricMV->GetTransform()->GetParameters()<<endl;
////    //metricMV->SetTransformParameters(parameters);
////    //metricMV->GetValue(parameters);
////    
////    
////    
////    OptimizerType::ScalesType scales(transform->GetNumberOfParameters());
////    //setting the scales and ranges for the parameters of the transform
////    const double translationScale = 2000.0; //dynamic range for tsl
////    const double rotationScale = 50.0; //dynamic range for rotations
////    
////    scales[0] = 1/rotationScale;
////    scales[1] = 1/rotationScale;
////    scales[2] = 1/rotationScale;
////    scales[3] = 1/translationScale;
////    scales[4] = 1/translationScale;
////    scales[5] = 1/translationScale;
////    
////    unsigned long numberOfIterations = 200;
////    //double gradientTolerance = 1e-4;
////    double valueTolerance = 1e-3;
////    double epsilonFunction = 1e-4;
////    
////    optimizer->SetCostFunction(metricMV);
////    optimizer->SetScales(scales);
////    optimizer->SetNumberOfIterations(numberOfIterations);
////    optimizer->SetValueTolerance(valueTolerance);
////    //optimizer->SetGradientTolerance(gradientTolerance);
////    optimizer->SetEpsilonFunction(epsilonFunction);
////    optimizer->SetInitialPosition(parameters);
////    
////    //optimizer->UseCostFunctionGradientOff();
////    
////    //observateur
////    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
////    optimizer->AddObserver(itk::IterationEvent(), observer);
////    
////    cout<<"start optimisation "<<endl;
////    //optimizer->StartOptimization();
////    
////    try {
////        optimizer->StartOptimization();
////    } catch (itk::ExceptionObject &e) {
////        cerr<<"error while optimizing cost function"<<endl;
////        cerr<<e<<endl;
////    }
////    
////    
////    cout<<"transform best param according to optimizer : "<<optimizer->GetCurrentPosition()<<endl;
////    cout<<"best param test transform : "<<transform->GetParameters()<<endl;
////    cout<<"best value : "<<optimizer->GetValue()<<endl;
////
//
//  //Amoeba optimizer
//    
//    registration->SetMetric(metric);
//    
//    
//    //setting des images
//    cout<<"setting images"<<endl;
//    registration->SetFixedImage(rescaled_IRM);
//    registration->SetMovingImage(US_shrunk);
//    //just to be sure for the metric
//    metric->SetFixed(rescaled_IRM);
//    metric->ComputeGradImage();
//    //metric->ComputeVesselnessImage();
//    metric->SetMoving(US_shrunk);
//    metric->ComputeMask();
//    
//    registration->SetFixedImageRegion(rescaled_IRM->GetBufferedRegion());
//    
//    
//    RegistrationType::ParametersType initialParameters = transform->GetParameters();
//    //setting des parametres optimizable
//    //test param BOBYQA porc1[-0.03514169419515513, 0.04780031905516767, 0.08510374453376125, -3.16002229212777, -0.7603111683282886, -1.5771778577960676]
//    initialParameters[0] = -0.03514169419515513;
//    initialParameters[1] = 0.04780031905516767;
//    initialParameters[2] = 0.08510374453376125;
//    initialParameters[3] = -3.16002229212777; //euler tsf = 6 parameters
//    initialParameters[4] = -0.7603111683282886;
//    initialParameters[5] = -1.5771778577960676;
//    
//    registration->SetInitialTransformParameters(initialParameters);
//    
//    cout<<"Initial transform param = "<<registration->GetTransform()->GetParameters()<<endl;
//    
//    /************
//     * OPTIMIZER
//     **************/
//    
//    //setting des params de l'optimizer
//    const unsigned int numberOfParameters = transform->GetNumberOfParameters();
//    AmoebaOptimizerType::ParametersType simplexDelta(numberOfParameters);
//    simplexDelta[0] =0.5;
//    simplexDelta[1] =0.5;
//    simplexDelta[2] =0.5;
//    simplexDelta[3] = 10;
//    simplexDelta[4] = 10;
//    simplexDelta[5] = 10;
//    
//    cout<<"verification simplex delta structure : "<<endl;
//    cout<<simplexDelta<<endl;
//    
//    optimizer->AutomaticInitialSimplexOff();
//    optimizer->SetInitialSimplexDelta(simplexDelta);
//    
//    //tolerance setting on param units and LC2 units
//    optimizer->SetParametersConvergenceTolerance(0.1);
//    optimizer->SetFunctionConvergenceTolerance(0.001);
//    
//    //max number of iterations to stop if no convergence
//    optimizer->SetMaximumNumberOfIterations(200);
//    
//    //on le met en maximisation
//    optimizer->SetMaximize(true);
//    
//    //Command observer
//    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
//    optimizer->AddObserver(itk::IterationEvent(), observer);
//    
//
//    cout<<"registration"<<endl;
//    
//        try {
//            registration->Update();
//            cout<<"Optimizer stop condition : "<<registration->GetOptimizer()->GetStopConditionDescription()<<endl;
//        } catch (itk::ExceptionObject &e) {
//            cerr<<"erreur dans le recalage"<<endl;
//            cerr<<e<<endl;
//            return EXIT_FAILURE;
//        }
//    
//
//    RegistrationType::ParametersType finalParameters = registration->GetLastTransformParameters();
//    
//    double bestValue = optimizer->GetValue();
//    
//    //Print out results
//    
//    cout<<"Results : "<<endl;
//    //cout<<"Translation vector : "<<"[ "<<TX<<", "<<TY<<", "<<TZ<<" ]"<<endl;
//    cout<<"parametres : "<<finalParameters<<endl;
//    cout<<"metric value : "<<bestValue<<endl;
//    
////    cout<<"formation de l'image recalee"<<endl;
////    
////    //TranslationType::Pointer finalTsl = TranslationType::New();
////    EulerTransformType::Pointer finalTsf = EulerTransformType::New();
////    
////    finalTsf->SetParameters(finalParameters);
////    finalTsf->SetFixedParameters(registration->GetTransform()->GetFixedParameters());
////    
////    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
////    resampler->SetTransform(finalTsf);
////    resampler->SetInput(image_US);
////    resampler->SetOutputDirection(finalTsf->GetInverseMatrix()*image_US->GetDirection());
////    resampler->SetSize(image_US->GetLargestPossibleRegion().GetSize());
////    resampler->SetOutputSpacing(image_US->GetSpacing());
////    resampler->SetOutputOrigin(finalTsf->GetInverseTransform()->TransformPoint(image_US->GetOrigin()));
////  
////    ImageType::Pointer outputImage = ImageType::New();
////    outputImage = resampler->GetOutput();
////    
////    WriterType::Pointer writer7 = WriterType::New();
////    string out7 =outputPath+"/finalregistreredUSAmoeba.nii.gz";
////    writer7->SetImageIO(io);
////    writer7->SetInput(outputImage);
////    writer7->SetFileName(out7);
////    try {
////        writer7->Update();
////    } catch (itk::ExceptionObject &e) {
////        cerr<<"error whilte writing registered image"<<endl;
////        cerr<<e<<endl;
////        return EXIT_FAILURE;
////    }
//    
    
    /***********************
     * TEST TSF DEFORMABLE
     **********************/
    
    BSplineTransformType::Pointer Btransform = BSplineTransformType::New();
    
    
    
    
    
    std::cout << "Done running routine!\n";
  
    tend = std::time(NULL);
    texec = std::difftime(tend,tbegin);
    
    cout<<"temps en s : "<<texec<<endl;
    return 0;
}
