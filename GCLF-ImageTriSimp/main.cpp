#include <QImage>
#include <QImageReader>
#include <QImageWriter>
#include "ParamParser.h"
#include "LinearSimplification/LinearSimplifier.h"
#include "CurvedSimplification/CurvedSimplifier.h"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string.hpp"

using namespace GCLF;
using namespace GCLF::Utils;
using namespace GCLF::Geometry;
using namespace GCLF::ImageTriSimp;
namespace bf = boost::filesystem;
namespace bj = boost::json;

void QImage2ImageT(const QImage& qimg, ImageT& img);
void ImageT2QImage(const ImageT& img, QImage& qimg);

int main(int argc, char* argv[])
{
  // args tips
  if (argc != 7)
  {
    printf("Only %d args provided.\n\n", argc);
    printf("Need args:\n");
    printf("arg[1]: parameters.\n");
    printf("input \"default\" to set default parameters or a json file to set parameters.\n");
    printf("arg[2]: input image path.\n");
    printf("arg[3]: input mesh path.\n");
    printf("arg[4]: output dir path.\n");
    printf("arg[5]: error bound.\n");
    printf("arg[6]: run mode.\n");
    return 1;
  }
  int param_idx = 1;
  Logger::InitLogger();

  // parse parameters
  std::string arg_param(argv[param_idx]);
  ParamTriangulator param;

  if (arg_param != "default")
  {
    bf::path json_file_path(argv[param_idx]);
    if (bf::is_regular_file(json_file_path))
    {
      fstream json_file;
      json_file.open(json_file_path.string(), fstream::in);
      if (json_file.is_open())
      {
        std::string json_str((std::istreambuf_iterator<char>(json_file)), std::istreambuf_iterator<char>());
        bj::stream_parser sp;
        sp.write(json_str.c_str());
        param.deserialize(sp.release().as_object());
        json_file.close();
      }
      else
      {
        Logger::user_logger->error("fail to open json file.");
        return 1;
      }
    }
    else
    {
      Logger::user_logger->error("wrong json file.");
      return 1;
    }
  }

  param_idx++;

  // parse input image
  bf::path in_image_path(argv[param_idx]);
  if (!bf::is_regular_file(in_image_path))
  {
    Logger::user_logger->error("error in input image.");
    return 1;
  }
  param_idx++;

  // parse input/output file/directory.
  bf::path in_model_path(argv[param_idx]);  param_idx++;
  bf::path out_data_path(argv[param_idx]);  param_idx++;

  if (!bf::is_regular_file(in_model_path))
  {
    Logger::user_logger->error("error in input model.");
    return 1;
  }

  if (!bf::is_directory(out_data_path))
  {
    Logger::user_logger->error("error output directory.");
    return 1;
  }

  if (arg_param == "default")
  {
    auto json_file_path = out_data_path;
    json_file_path.append("config.json");
    if (!bf::is_regular_file(json_file_path))
    {
      fstream json_file;
      json_file.open(json_file_path.string(), fstream::out);
      if (json_file.is_open())
      {
        auto json_obj = param.serialize();
        auto json_str = bj::serialize(json_obj);
        json_file << json_str;
        json_file.close();
      }
    }
  }

  // parse error bound
  double error_bound;
  try
  {
    std::string error_str(argv[param_idx]);
    param_idx++;
    error_bound = std::stod(error_str);
  }
  catch (...)
  {
    Logger::user_logger->error("error in parsing error bound.");
    return 1;
  }
  param.error_bound = error_bound;

  // parse target vertices numbers.
  std::string mode(argv[param_idx]);
  param_idx++;

  // read input image
  ImageT input_image;
  QImage qt_input_image;
  QImageReader qt_image_reader;
  QString qt_image_path = QString::fromStdString(in_image_path.string());

  qt_image_reader.setFileName(qt_image_path);
  qt_image_reader.setAutoDetectImageFormat(true);
  if (!qt_image_reader.canRead())
  {
    Logger::user_logger->error("error {} in reading input image.", qt_image_reader.error());
    return 1;
  }
  qt_image_reader.read(&qt_input_image);
  QImage2ImageT(qt_input_image, input_image);

  // read input mesh
  LMeshT in_linear_mesh;  bool get_linear = false;
  BMeshT in_curved_mesh;  bool get_curved = false;

  if (in_model_path.extension().string() == ".obj")
  {
    OpenMesh::IO::read_mesh(in_linear_mesh, in_model_path.string());
    get_linear = true;
  }
  else if (in_model_path.extension().string() == ".bma")
  {
    in_curved_mesh.read_ascii(in_model_path.string());
    get_curved = true;
  }
  else
  {
    Logger::user_logger->error("unrecognized model type.");
  }

  if (mode == "split")
  {
    // set logger file
    auto log_file_path = out_data_path;
    log_file_path.append("split_log.txt");
    Logger::updateFileLog(true, spdlog::level::info, log_file_path.string());

    LinearSimplifier simplifier;

    simplifier.initialize(input_image, in_linear_mesh, param);
    simplifier.set_error_bound(param.error_bound);

    // save out mesh
    auto out_mesh_path = out_data_path;
    out_mesh_path.append("bounded_linear_mesh.obj");
    OpenMesh::IO::write_mesh(*simplifier.mesh, out_mesh_path.string());
    // save out image
    auto out_image_path = out_data_path;
    out_image_path.append("bounded_linear_image.png");
    QImage qt_out_image = qt_input_image;
    ImageT2QImage(*simplifier.out_image, qt_out_image);
    QImageWriter imageWriter;
    imageWriter.setFileName(QString::fromStdString(out_image_path.string()));
    imageWriter.write(qt_out_image);
  }
  else if (mode == "linear")
  {
    // set logger file
    auto log_file_path = out_data_path;
    log_file_path.append("linear_log.txt");
    Logger::updateFileLog(true, spdlog::level::info, log_file_path.string());

    LinearSimplifier simplifier;

    simplifier.initialize(input_image, in_linear_mesh, param);
    simplifier.set_error_bound(param.error_bound);
    simplifier.simplify();

    // save out mesh
    auto out_mesh_path = out_data_path;
    out_mesh_path.append("simplified_linear_mesh.obj");
    OpenMesh::IO::write_mesh(*simplifier.mesh, out_mesh_path.string());
    // save out image
    auto out_image_path = out_data_path;
    out_image_path.append("simplified_linear_image.png");
    QImage qt_out_image = qt_input_image;
    ImageT2QImage(*simplifier.out_image, qt_out_image);
    QImageWriter imageWriter;
    imageWriter.setFileName(QString::fromStdString(out_image_path.string()));
    imageWriter.write(qt_out_image);
  }
  else if (mode == "curved")
  {
    // set logger file
    auto log_file_path = out_data_path;
    log_file_path.append("curved_log.txt");
    Logger::updateFileLog(true, spdlog::level::info, log_file_path.string());

    ASSERT(get_curved || get_linear, "fail to get mesh.");
    CurvedSimplifier simplifier;

    if (get_linear)
      simplifier.initialize(input_image, in_linear_mesh, param);
    else if (get_curved)
      simplifier.initialize(input_image, in_curved_mesh, param);
    simplifier.set_error_bound(param.error_bound);
    simplifier.simplify();

    // save out mesh
    auto out_mesh_path = out_data_path;
    out_mesh_path.append("simplified_curved_mesh.bma");
    simplifier.mesh->write_ascii(out_mesh_path.string());
    // save out image
    auto out_image_path = out_data_path;
    out_image_path.append("simplified_curved_image.png");
    QImage qt_out_image = qt_input_image;
    ImageT2QImage(*simplifier.out_image, qt_out_image);
    QImageWriter imageWriter;
    imageWriter.setFileName(QString::fromStdString(out_image_path.string()));
    imageWriter.write(qt_out_image);
  }
  else
  {
    Logger::user_logger->error("unimplemented mode.");
    return 1;
  }

  return 0;
}


void QImage2ImageT(const QImage& qimg, ImageT& img)
{
  img.resize(qimg.width(), qimg.height());
  // QImage starts from Upper left.
  for (int y = 0;y < img.height;y++)
    for (int x = 0;x < img.width;x++)
    {
      QRgb qcolor = qimg.pixel(x, y);
      ImageT::Color color;
    #if defined(CHANNEL_0_1)
      ImageT::red_channel(color) = qRed(qcolor) / 255.;
      ImageT::green_channel(color) = qGreen(qcolor) / 255.;
      ImageT::blue_channel(color) = qBlue(qcolor) / 255.;
    #elif defined(CHANNEL_0_255)
      ImageT::red_channel(color) = qRed(qcolor);
      ImageT::green_channel(color) = qGreen(qcolor);
      ImageT::blue_channel(color) = qBlue(qcolor);
    #endif
      img.pixel_color(x, y) = color;
    }
}

void ImageT2QImage(const ImageT& img, QImage& qimg)
{
  qimg.scaled(img.width, img.height);
  for (int y = 0;y < img.height;y++)
    for (int x = 0;x < img.width;x++)
    {
      QRgb qcolor = 0;
      ImageT::Color color = img.pixel_color(x, y);
    #if defined(CHANNEL_0_1)
      uint32_t qred = (int)round(ImageT::red_channel(color) * 255.f);
      uint32_t qgreen = (int)round(ImageT::green_channel(color) * 255.f);
      uint32_t qblue = (int)round(ImageT::blue_channel(color) * 255.f);
    #elif defined(CHANNEL_0_255)
      uint32_t qred = (int)round(ImageT::red_channel(color));
      uint32_t qgreen = (int)round(ImageT::green_channel(color));
      uint32_t qblue = (int)round(ImageT::blue_channel(color));
    #endif
      qcolor = (255u << 24) | (qred << 16) | (qgreen << 8) | qblue;
      qimg.setPixel(x, y, qcolor);
    }
}