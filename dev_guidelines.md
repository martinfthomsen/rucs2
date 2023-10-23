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
1. Create a new development branch from the master branch
   a) if it is a new patch, prefix it with the release type: "patch/"
   b) if it is a new feature, prefix it with the release type: "feature/"
   c) if it is a new product, prefix it with the release type: "product/"
2. Name the release using a short but descriptive title in snake case.
fx feature/add_explore_mode
3. Make the needed changes to your local the repo
4. Build the Docker image
   docker-compose build
5. Test and verify that the code functions as expected
6. Create pull request for the branch to the master
    - Go to the pull request pane for the repo https://github.com/martinfthomsen/rucs2/pulls
    - Click "New pull request" button
    - Select the development branch you wish to merge with the master branch
    - Verify everything look ok for merging and resolve any listed issues
    - Click "Create pull request"
    - Write some description of the new additions/changes
    - Make sure to link to the issues which this branch resolves (using the # tag fx. #1 for issue nummer one 1)
    - Submit the pull request
7. Merge pull request
    - Go to the pull request pane for the repo https://github.com/martinfthomsen/rucs2/pulls
    - Click on the pull request you wish to merge
    - Verify that everything is in order for the merger (or reject it)
    - Click to merge the pull request
    - Confirm the merger
8. Publish a new release for the new commit to the master branch
    - Go to the newly merged commit to the master branch (fx on your local instance)
    - Add tag with a name according to the release using semantic versioning. fx. v1.2.0. (remember to push it to the server if you add the tag locally)
      The first number is the product release
      The second number is the feature release
      The third number is the patch release
      According to the branch type being merged, increment the correct number.
      If it is a new feature release, remember to reset the patch number to 0.
      And if it is a product release, both feature and patch number should be reset to 0.
    - Go to the releases pane for the repo https://github.com/martinfthomsen/rucs2/releases
    - Click on "Draft a new release"
    - Choose the tag of the new release
    - Use the tag name as the release name
    - Click "Generate release notes"
    - Finish the release notes with ample documentation of the changes
      Suggestion: Look up the pull request description, as they might contain a good description of the changes.
      Suggestion: Look at a previos release note for inspiration of how to write the new release note
      Suggestion: Look at the contained git commits descriptions for further details
    - Publish release


## Pushing to DockerHub ##
1. login to docker from cmdline:
```
docker login --username=<USER>
# provide password
```
2. Checkout the git tag you wish to submit
```
git checkout tags/<TAG>
version=$(git tag -l --points-at HEAD)
```
3. Build image
```
docker-compose build
```
4. Add local docker tags for the build
```
docker tag rucs <USER>/rucs2:$version
docker tag rucs <USER>/rucs2
```
5. Push tagged builds to DockerHub
```
docker push <USER>/rucs2:$version
docker push <USER>/rucs2
```

### Example ###
```
docker login --username=mcft
git checkout tags/v1.3.0
version=$(git tag -l --points-at HEAD)
docker-compose build
docker tag rucs mcft/rucs2:$version
docker tag rucs mcft/rucs2
docker push mcft/rucs2:$version
docker push mcft/rucs2
```


## Debugging ##
The easiest way of debugging the code is to enter the container in interactive
mode and run the code interactively, with quiet and clean_run options turned off.

### Steps ###
1. Create test directory
2. Copy your test data into your test directory (positives and negatives subdirs)
   (or download them)
```
#!bash
docker run -it --entrypoint download_genomes.sh --rm -v `pwd`/positives:/workdir rucs CP000672.1 ARBW00000000
docker run -it --entrypoint download_genomes.sh --rm -v `pwd`/negatives:/workdir rucs JWIZ01
```
3. Open the RUCS container and initiate
```
#!bash
docker run -it --entrypoint python3 --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb rucs
import sys
sys.path.append('/tools/')
from primer_core_tools import *
load_global_settings('settings.default.cjson')
positives = [x for x in glob.glob("positives/*") if check_file_type(x) == 'fasta']
negatives = [x for x in glob.glob("negatives/*") if check_file_type(x) == 'fasta']
reference = positives[0]
```
4. Run the program you wish to test/debug
#### Full run: ####
```
#!bash
main(positives, negatives, reference, quiet=False, clean_run=False, annotate=True)
```
#### Explore run: ####
```
#!bash
explore(positives, negatives, quiet=False, clean_run=False, reuse=True)
```

### Access the container through bash ###
Sometimes what you need to test lays outside the primer_core_tools python library.
For instance in the case of installing new version of the blast tools or similar.
For this, run the following:
```
#!bash
docker run -it --entrypoint bash --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb rucs
```
Subsequent runs of RUCS can be achieved by running python3 and running the above
commands from step 3 and 4.
