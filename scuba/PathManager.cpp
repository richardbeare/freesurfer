/**
 * @file  PathManager.cpp
 * @brief Scuba manager for Path<float>s
 *
 * Managers multiple Path<floats>, relaying pathChanged messages from
 * them to its Listener. Handles Tcl commands to read and write path
 * files.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/16 20:18:30 $
 *    $Revision: 1.12 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <fstream>
#include "PathManager.h"

using namespace std;

PathManager&
PathManager::GetManager() {

  static PathManager* sManager = NULL;
  if ( NULL == sManager ) {
    sManager = new PathManager();

    TclCommandManager& commandMgr = TclCommandManager::GetManager();
    commandMgr.AddCommand( *sManager, "WritePathFile", 1, "fileName",
                           "Write paths to a file." );
    commandMgr.AddCommand( *sManager, "ReadPathFile", 1, "fileName",
                           "Read paths from a file." );
  }

  return *sManager;
}

PathManager::PathManager () :
    Broadcaster( "PathManager" ),
    Listener( "PathManager" ) {
}

void
PathManager::ManagePath ( Path<float>& iPath ) {

  iPath.AddListener( *this );
  mPaths.push_back( &iPath );

  // Notify listeners of changee.
  int pathID = iPath.GetID();
  SendBroadcast( "pathChanged", (void*)&pathID );
}

void
PathManager::UnmanagePath ( Path<float> const& iPath ) {

  // Search for a path with the same ID and remove it from our list.
  vector<Path<float>*>::iterator tPath;
  for ( tPath = mPaths.begin(); tPath != mPaths.end(); ++tPath ) {
    Path<float>* path = *tPath;
    if ( path->GetID() == iPath.GetID() ) {
      int pathID = path->GetID();
      mPaths.erase( tPath );

      // Notify listeners of changee.
      SendBroadcast( "pathChanged", (void*)&pathID );
      break;
    }
  }
}

vector<Path<float>*> const&
PathManager::GetPathList () const {
  return mPaths;
}

vector<Path<float>*>&
PathManager::GetPathList () {
  return mPaths;
}

TclCommandManager::TclCommandResult
PathManager::DoListenToTclCommand ( char* isCommand, int,
                                    char** iasArgv ) {

  // WritePathFile
  if ( 0 == strcmp( isCommand, "WritePathFile" ) ) {

    string fnPaths = iasArgv[1];
    WritePathFile( fnPaths );
  }

  // ReadPathFile
  if ( 0 == strcmp( isCommand, "ReadPathFile" ) ) {

    string fnPaths = iasArgv[1];
    ReadPathFile( fnPaths );
  }


  return ok;
}

void
PathManager::DoListenToMessage ( string isMessage, void* iData ) {

  if ( isMessage == "pathChanged" ||
       isMessage == "pathVertexAdded" ) {

      SendBroadcast( isMessage, iData );
  }
}

void
PathManager::ReadPathFile ( string const& ifnPaths ) {

  // Try to open the file.
  ifstream fPaths( ifnPaths.c_str(), ios::in );
  if ( !fPaths || fPaths.bad() ) {
    throw runtime_error( "Couldn't open paths file." );
  }

  // Read the number of paths.
  int cPaths;
  fPaths >> cPaths;

  // For each one, create a new path and have it read from the
  // stream. Then manage that path.
  for ( int nPath = 0; nPath < cPaths; nPath++ ) {

    Path<float>* path = new Path<float>;
    path->ReadFromStream( fPaths );
    ManagePath( *path );

    // Notify listeners of changee.
    int pathID = path->GetID();
    SendBroadcast( "pathChanged", (void*)&pathID );
  }

  fPaths.close();
}

void
PathManager::WritePathFile ( string const& ifnPaths ) {

  // Open the file.
  ofstream fPaths( ifnPaths.c_str(), ios::out );
  if ( !fPaths || fPaths.bad() ) {
    throw runtime_error( "Couldn't write paths file." );
  }

  // Write the number of paths.
  fPaths << mPaths.size() << endl;

  // For each of our paths, have the path write to the stream.
  vector<Path<float>*>::iterator tPath;
  for ( tPath = mPaths.begin(); tPath != mPaths.end(); ++tPath ) {
    Path<float>* path = *tPath;
    path->WriteToStream( fPaths );
  }

  fPaths.close();
}
