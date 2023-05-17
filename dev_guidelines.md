# Procedures and guidelines for developer #
Based on GitHub Flow.

All commits to the master branch are new releases. There are three release types
1. Product releases, are releases that break backwards compatibility or which are
so huge they require a major release
2. Feature releases, are releases that are backwards compatible, such as addition
of new features and functions, etc. These can often contain bug fixes and other
minor changes to the repo content.
3. Patch releases, are a collection of hot-fixes that fix a bug in the code,
refactoring of code that does not impact the functionality, aesthetic changes,
and text updates and corrections to the readme, etc. which cannot wait for the
next minor code release.

## Making a Release ##
To make a new release:
1. Create a new branch from the master branch
   a) if it is a new patch, prefix it with the release type: "patch/"
   b) if it is a new feature, prefix it with the release type: "feature/"
   c) if it is a new product, prefix it with the release type: "product/"
2. Name the release using a short but descriptive title in snake-case.
fx feature/add_explore_mode
3. Make the needed changes to the codebase
4. Build the Docker image and test the code functions as expected
5. Create pull request for the branch to the master
6. Merge pull request
7. Publish a new release for the new commit to the master branch
