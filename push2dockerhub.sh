DOCKER HUB

create user (online - https://hub.docker.com/)
create repository (online - https://hub.docker.com/)
docker tag <IMAGE> <USER>/<REPOSITORY>:<VERSION>
docker push <USER>/<REPOSITORY>:<VERSION>


# Get current version from git tag of current commit
git checkout tags/1.0.1
version=$(git tag -l --points-at HEAD)
pv=${version%.*}
ppv=${pv%.*}
docker-compose build

# Update all versions variants
docker tag rucs genomicepidemiology/rucs:$version
docker tag rucs genomicepidemiology/rucs:$pv
docker tag rucs genomicepidemiology/rucs:$ppv
docker tag rucs genomicepidemiology/rucs

# Push all versions variants
docker push genomicepidemiology/rucs:$version
docker push genomicepidemiology/rucs:$pv
docker push genomicepidemiology/rucs:$ppv
docker push genomicepidemiology/rucs


HotFix procedure
==================
# Create Hotfix branch (increment previous hotfix version)
# Make changes to the files
# Build and test the docker image
docker-compose build
docker tag genomicepidemiology/rucs rucs
docker run --rm -v `pwd`:/workdir -v $BLASTDB:/blastdb rucs test
# Stage and Commit changes, when everything is in order
# Finalise Hotfix branch
# Checkout the version tag
# Update Docker tags
	# Get current version from git tag of current commit
	git checkout tags/1.0.2
	version=$(git tag -l --points-at HEAD)
	pv=${version%.*}
	ppv=${pv%.*}
	# Update all versions tags
docker tag rucs genomicepidemiology/rucs:$version
docker tag rucs genomicepidemiology/rucs:$pv
docker tag rucs genomicepidemiology/rucs:$ppv
docker tag rucs genomicepidemiology/rucs

# Push Docker tags
docker push genomicepidemiology/rucs:$version
docker push genomicepidemiology/rucs:$pv
docker push genomicepidemiology/rucs:$ppv
docker push genomicepidemiology/rucs

# Clean out unused images



Made changes to wrong branch, and want to switch branch without loosing the changes
=======================
# Stash your changes (add optional message)
# Create or which to new branch
# Apply stash to the new branch
# Delete hunks which you don't want to apply to the new branch
# Stage and commit your changes
# Delete stash when you are finished with it



