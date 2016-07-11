#include <gtest/gtest.h>
#include <easylogging++.h>
#include "./test_models.cpp"
#include "./test_loadvcf.cpp"
#include "./test_liabilitymodel.cpp"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
