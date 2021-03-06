set(TestXML_ARGS "DATA{${VTK_TEST_INPUT_DIR}/sample.xml}")
vtk_add_test_cxx(vtkIOXMLCxxTests tests
  TestAMRXMLIO.cxx,NO_VALID
  TestDataObjectXMLIO.cxx,NO_VALID
  TestReadDuplicateDataArrayNames.cxx,NO_DATA,NO_VALID
  TestXML.cxx,NO_DATA,NO_VALID,NO_OUTPUT
  TestXMLGhostCellsImport.cxx
  TestXMLHierarchicalBoxDataFileConverter.cxx,NO_VALID
  TestXMLHyperTreeGridIO.cxx,NO_VALID
  TestXMLMappedUnstructuredGridIO.cxx,NO_DATA,NO_VALID
  TestXMLToString.cxx,NO_DATA,NO_VALID,NO_OUTPUT
  TestXMLUnstructuredGridReader.cxx
  TestXMLWriterWithDataArrayFallback.cxx,NO_VALID
  )

# Each of these most be added in a separate vtk_add_test_cxx
vtk_add_test_cxx(vtkIOXMLCxxTests tests
  TestXMLCompositeDataReaderDistribution.cxx,NO_VALID,NO_OUTPUT
  "DATA{${VTK_TEST_INPUT_DIR}/distTest.vtm,REGEX:distTest_[0-9]_0.vtp}"
  )
vtk_add_test_cxx(vtkIOXMLCxxTests tests
  TestXMLReaderBadImageData,TestXMLReaderBadData.cxx,NO_VALID,NO_OUTPUT "DATA{${VTK_TEST_INPUT_DIR}/badImageData.xml}"
)
vtk_add_test_cxx(vtkIOXMLCxxTests tests
  TestXMLReaderBadPolyData,TestXMLReaderBadData.cxx,NO_VALID,NO_OUTPUT "DATA{${VTK_TEST_INPUT_DIR}/badPolyData.xml}"
)
vtk_add_test_cxx(vtkIOXMLCxxTests tests
  TestXMLReaderBadRectilinearGridData,TestXMLReaderBadData.cxx,NO_VALID,NO_OUTPUT "DATA{${VTK_TEST_INPUT_DIR}/badRectilinearGridData.xml}"
)
vtk_add_test_cxx(vtkIOXMLCxxTests tests
  TestXMLReaderBadUnstucturedGridData,TestXMLReaderBadData.cxx,NO_VALID,NO_OUTPUT "DATA{${VTK_TEST_INPUT_DIR}/badUnstructuredGridData.xml}"
)
vtk_add_test_cxx(vtkIOXMLCxxTests tests
  TestXMLReaderBadUniformGridData,TestXMLReaderBadData.cxx,NO_VALID,NO_OUTPUT "DATA{${VTK_TEST_INPUT_DIR}/badUniformGridData.xml}"
)

vtk_test_cxx_executable(vtkIOXMLCxxTests tests)

add_executable(TestXMLCInterface MACOSX_BUNDLE TestXMLCInterface.c)
target_link_libraries(TestXMLCInterface LINK_PRIVATE vtkIOXML)
add_test(NAME TestXMLCInterface COMMAND TestXMLCInterface)
