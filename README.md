# Cinch Example Project

Congratulations! If you are reading this file, you have successfully created
a cinch skeleton project.  From here, you should try to configure and build
the example project:

    % mkdir build
    % cd build
    % cmake -DENABLE\_DOCUMENTATION=ON -DENABLE\_DOXYGEN=ON -DENABLE\_JENKINS\_OUTPUT -DENABLE\_UNIT\_TESTS
    % make

After you have built the project, you should have several executable and
documentation files in your build tree:

    bin/app
    lib/libexample.a
    test/example/sanity
    doc/user-guide-0.0.pdf
    doc/developer-guide-0.0.pdf
    doc/doxygen

To run the application example:

    % bin/app

This should produce, "Hello World".

To run the unit tests:

    % make test

Or optionally:

    % ctest -N (list the available unit tests)
    % ctest -V -R sanity (run tests verbose that match 'sanity')

# Push your new project to a remote git repository

To push to a remote, you will first need to create the repository on
whichever git server you plan to use.  Generally, this requires
logging-in to the server and creating an empty project.  The server
will usually tell you the path to the new project, which might look
something like:

    git@github.com:organization/project.git

Using this url as an example, to push your cinch project to the remote:

    % cd /path/to/project
    % git remote add origin git@github.com:organization/project.git
    % git push -u origin master